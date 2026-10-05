--[[
response-letter.lua

One semantic model, N renderers.

Semantic layer (format-agnostic, owned by this filter):

  letterhead → sections → response blocks

  ::: {.comment}            reviewer/editor comment, auto-numbered
  ::: {.comment .general}   unnumbered comment ("overall assessment")
  ::: {.reply}              authors' reply; optional status="done|partial|declined|todo"
  ::: {.changes}            list of manuscript changes
  ::: {.quoted loc="..."}   verbatim quote from the revised manuscript
  ::: {.response-summary}   replaced by the collected summary of changes

Sectioning: a plain level-1 header starts the next reviewer (counter R);
`# ... {.editor}` starts the editor section (counter E); `{.unnumbered}`
headers leave the counters alone. Comments are numbered R.C / E.N and get
the anchor #cmt-R-C / #cmt-E-N, so "[Comment 1.2](#cmt-1-2)" links work in
every format.

Renderers (adding a format = adding one entry in `render`):
  latex → tcolorbox environments from header.tex; the color palette and the
          fancyhdr running header are injected into the preamble from metadata
  docx  → custom paragraph styles from reference.docx; the running header
          lives in reference.docx (STYLEREF/TITLE fields)
  html  → CSS classes; the palette is injected as CSS custom properties

Letter metadata is read from the `response-letter` YAML key (see README);
the letterhead and closing are built as Pandoc AST blocks so a single code
path serves all formats.
]]

local stringify = pandoc.utils.stringify

-------------------------------------------------------------------- config

local KINDS = { "comment", "reply", "changes", "quoted" }

-- Defensive fallbacks only; canonical defaults live in _extension.yml.
local cfg = {
  labels = {
    comment  = "Comment",
    reply    = "Reply",
    changes  = "Changes",
    quoted   = "Manuscript",
    reviewer = "Reviewer",
    editor   = "Editor",
    summary  = "Summary of Changes",
    re       = "Re:",
    dear     = "Dear",
    round    = "Revision",
  },
  status = {
    done     = "✓",
    partial  = "partially addressed",
    declined = "respectfully declined",
    todo     = "TO DO",
  },
  colors = {
    comment = "#B22222",
    reply   = "#1E90FF",
    changes = "#8FBC8F",
    quoted  = "#696969",
  },
  -- HTML background tints: fraction of the base color mixed into white,
  -- chosen to reproduce the v0.1.0 palette.
  tints = { comment = 0.05, reply = 0.04, changes = 0.12, quoted = 0.07 },
  letter = {},
  signature = nil,          -- raw meta value (inlines or blocks)
  shorttitle = nil,
  running_header = true,
}

local LANG_KEYS = {
  comment  = "response-comment-label",
  reply    = "response-reply-label",
  changes  = "response-changes-label",
  quoted   = "response-quoted-label",
  reviewer = "response-reviewer-label",
  editor   = "response-editor-label",
  summary  = "response-summary-title",
  re       = "response-re-label",
  dear     = "response-dear-label",
  round    = "response-round-label",
}

local STATUS_KEYS = {
  done     = "response-status-done",
  partial  = "response-status-partial",
  declined = "response-status-declined",
  todo     = "response-status-todo",
}

local DOCX_LABEL_STYLE = "Response Label"

local latex_env = {
  comment = "CommentBox",
  reply   = "ReplyBox",
  changes = "ChangesBox",
  quoted  = "QuoteBox",
}

local docx_style = {
  comment = "Response Comment",
  reply   = "Response Reply",
  changes = "Response Changes",
  quoted  = "Response Quote",
}

local typst_fn = {
  comment = "rl-comment",
  reply   = "rl-reply",
  changes = "rl-changes",
  quoted  = "rl-quoted",
}

local function warn(msg)
  if quarto and quarto.log and quarto.log.warning then
    quarto.log.warning("response-letter: " .. msg)
  else
    io.stderr:write("[WARNING] response-letter: " .. msg .. "\n")
  end
end

local function meta_str(v)
  if v == nil then return nil end
  return stringify(v)
end

local function read_meta(meta)
  local lang = meta.language
  if lang then
    for key, langkey in pairs(LANG_KEYS) do
      if lang[langkey] then cfg.labels[key] = stringify(lang[langkey]) end
    end
    for key, langkey in pairs(STATUS_KEYS) do
      if lang[langkey] then cfg.status[key] = stringify(lang[langkey]) end
    end
  end

  local rl = meta["response-letter"]
  if rl then
    local L = cfg.letter
    L.journal = meta_str(rl.journal)
    if rl.manuscript then
      L.id       = meta_str(rl.manuscript.id)
      L.mtitle   = meta_str(rl.manuscript.title)
      L.mversion = meta_str(rl.manuscript.version)
    end
    L.round = meta_str(rl.round)
    if rl.editor then
      L.editor_name = meta_str(rl.editor.name)
      L.editor_role = meta_str(rl.editor.role)
    end
    L.salutation = meta_str(rl.salutation)
    L.closing    = meta_str(rl.closing)
    cfg.signature = rl.signature
    if rl.colors then
      for _, kind in ipairs(KINDS) do
        if rl.colors[kind] then cfg.colors[kind] = stringify(rl.colors[kind]) end
      end
    end
    if rl["running-header"] ~= nil then
      local v = rl["running-header"]
      cfg.running_header = (v == true) or (meta_str(v) == "true")
    end
    -- tri-state: true inserts the summary before the first section when no
    -- .response-summary div is present; false suppresses summaries even if
    -- a div exists; nil (unset) leaves placement to the div alone.
    if rl.summary ~= nil then
      local v = rl.summary
      if type(v) == "boolean" then
        L.summary = v
      else
        L.summary = (meta_str(v) == "true")
      end
    end
  end

  cfg.shorttitle = meta_str(meta.shorttitle) or meta_str(meta.title)
end

-------------------------------------------------------------- target format

local function target_format()
  if quarto and quarto.doc then
    if quarto.doc.is_format("latex") then return "latex" end
    if quarto.doc.is_format("docx") then return "docx" end
    if quarto.doc.is_format("typst") then return "typst" end
    return "html"
  end
  if FORMAT:match("latex") then return "latex" end
  if FORMAT:match("docx") then return "docx" end
  if FORMAT:match("typst") then return "typst" end
  return "html"
end

local target = target_format()

--------------------------------------------------------------------- state

local sec = { kind = nil }            -- nil | "reviewer" | "editor"
local reviewer_n = 0
local comment_n  = 0                  -- per-reviewer counter
local editor_n   = 0                  -- global editor-comment counter
local summary_items = pandoc.List()   -- { label, anchor, summary, status }
local todo_items    = pandoc.List()
local last_comment  = nil             -- summary item of the latest numbered comment
local warned_unsectioned = false

------------------------------------------------------------------- helpers

local function latex_escape(s)
  return (s:gsub("[%%&#_{}$^~\\]", {
    ["%"] = "\\%%", ["&"] = "\\&", ["#"] = "\\#", ["_"] = "\\_",
    ["{"] = "\\{", ["}"] = "\\}", ["$"] = "\\$",
    ["^"] = "\\textasciicircum{}", ["~"] = "\\textasciitilde{}",
    ["\\"] = "\\textbackslash{}",
  }))
end

-- Escape a string for a Typst double-quoted string literal.
local function typst_str(s)
  return (s:gsub("\\", "\\\\"):gsub('"', '\\"'))
end

-- "Reviewer 2" / "Editor": the text shown in the running header. LaTeX sets
-- \rightmark; Typst drops a hidden <rl-mark> the page header queries for the
-- current section (an empty text clears it, matching \markright{} for
-- .unnumbered sections).
local function section_mark(text)
  if not cfg.running_header then return nil end
  if target == "latex" then
    return pandoc.RawBlock("latex", "\\markright{" .. latex_escape(text) .. "}")
  elseif target == "typst" then
    return pandoc.RawBlock("typst",
      string.format('#metadata("%s")<rl-mark>', typst_str(text)))
  end
  return nil
end

local function comment_number()
  if sec.kind == "editor" then
    editor_n = editor_n + 1
    return string.format("E.%d", editor_n), string.format("cmt-E-%d", editor_n)
  end
  comment_n = comment_n + 1
  return string.format("%d.%d", reviewer_n, comment_n),
         string.format("cmt-%d-%d", reviewer_n, comment_n)
end

local function label_text(kind, div, number)
  if kind == "comment" then
    if number then
      return string.format("%s %s", cfg.labels.comment, number)
    end
    return cfg.labels.comment .. ":"
  elseif kind == "quoted" then
    local loc = div.attributes["loc"]
    if loc and loc ~= "" then
      return string.format("%s (%s):", cfg.labels.quoted, loc)
    end
    return cfg.labels.quoted .. ":"
  elseif kind == "reply" then
    local status = div.attributes["status"]
    if status == "done" then
      return string.format("%s %s:", cfg.labels.reply, cfg.status.done)
    elseif status == "partial" or status == "declined" then
      return string.format("%s (%s):", cfg.labels.reply, cfg.status[status])
    elseif status == "todo" then
      return string.format("%s (%s):", cfg.labels.reply, cfg.status.todo)
    elseif status then
      warn(string.format(
        "unknown reply status '%s' (expected done|partial|declined|todo)", status))
    end
    return cfg.labels.reply .. ":"
  end
  return cfg.labels[kind] .. ":"
end

-- Text → Inlines; in LaTeX output "✓" is replaced by pifont's \ding{51}
-- because the default text and math fonts silently drop the Unicode glyph
-- (unicode-math remaps \checkmark to a missing math-font slot under
-- lualatex/xelatex; Zapf Dingbats works under every engine).
local function text_inlines(text)
  if target ~= "latex" or not text:find("✓", 1, true) then
    return pandoc.Inlines({ pandoc.Str(text) })
  end
  local inlines = pandoc.Inlines({})
  local pos = 1
  while true do
    local s, e = text:find("✓", pos, true)
    if not s then
      if pos <= #text then inlines:insert(pandoc.Str(text:sub(pos))) end
      break
    end
    if s > pos then inlines:insert(pandoc.Str(text:sub(pos, s - 1))) end
    inlines:insert(pandoc.RawInline("latex", "\\ding{51}"))
    pos = e + 1
  end
  return inlines
end

local function label_block(kind, text, status)
  local para = pandoc.Para({ pandoc.Strong(text_inlines(text)) })
  if target == "docx" then
    return pandoc.Div({ para },
      pandoc.Attr("", {}, { ["custom-style"] = DOCX_LABEL_STYLE }))
  end
  local classes = { "box-label", "box-label-" .. kind }
  if status then table.insert(classes, "box-label-" .. status) end
  return pandoc.Div({ para }, pandoc.Attr("", classes))
end

----------------------------------------------------------------- renderers

-- Each renderer: (kind, div, content) -> pandoc.List of blocks.
local render = {}

render.latex = function(kind, div, content)
  local blocks = pandoc.List()
  if div.identifier ~= "" then
    -- Empty span carries the link target; pandoc emits the anchor syntax
    -- matching its own internal-link output.
    blocks:insert(pandoc.Plain({
      pandoc.Span({}, pandoc.Attr(div.identifier))
    }))
  end
  blocks:insert(pandoc.RawBlock("latex", "\\begin{" .. latex_env[kind] .. "}"))
  blocks:extend(content)
  blocks:insert(pandoc.RawBlock("latex", "\\end{" .. latex_env[kind] .. "}"))
  return blocks
end

render.typst = function(kind, div, content)
  local blocks = pandoc.List()
  blocks:insert(pandoc.RawBlock("typst", "#" .. typst_fn[kind] .. "["))
  blocks:extend(content)
  -- A Typst label attaches to the PRECEDING element, so the cross-ref anchor
  -- must follow the box. (The empty-span-before-the-box trick render.latex
  -- uses yields an unattached label and a fatal "label does not exist" error
  -- when something #link()s to it.)
  if div.identifier ~= "" then
    blocks:insert(pandoc.RawBlock("typst", "]\n<" .. div.identifier .. ">"))
  else
    blocks:insert(pandoc.RawBlock("typst", "]"))
  end
  return blocks
end

-- pandoc cannot carry paragraph styles into real Word list items (they
-- always get pandoc's own list style, losing the box shading), so lists
-- inside boxes become one paragraph per item with a literal bullet/number,
-- using the "<box style> List" styles defined in reference.docx (hanging
-- indent, inherited border and shading).
local function docx_listify(blocks, style, level)
  level = level or 0
  local indent = string.rep("\u{2003}", level)
  local out = pandoc.Blocks({})

  local function listed_item(prefix, item)
    local item_out = pandoc.Blocks({})
    local first = true
    for _, b in ipairs(item) do
      if b.t == "BulletList" or b.t == "OrderedList" then
        item_out:extend(docx_listify(pandoc.Blocks({ b }), style, level + 1))
      elseif first and (b.t == "Plain" or b.t == "Para") then
        first = false
        local inlines = pandoc.Inlines({ pandoc.Str(prefix) })
        inlines:extend(b.content)
        item_out:insert(pandoc.Div({ pandoc.Para(inlines) },
          pandoc.Attr("", {}, { ["custom-style"] = style })))
      else
        item_out:insert(b)
      end
    end
    return item_out
  end

  for _, block in ipairs(blocks) do
    if block.t == "BulletList" then
      for _, item in ipairs(block.content) do
        out:extend(listed_item(indent .. "• ", item))
      end
    elseif block.t == "OrderedList" then
      local n = block.start or 1
      for _, item in ipairs(block.content) do
        out:extend(listed_item(indent .. tostring(n) .. ". ", item))
        n = n + 1
      end
    else
      out:insert(block)
    end
  end
  return out
end

render.docx = function(kind, div, content)
  content = docx_listify(content, docx_style[kind] .. " List")
  local styled = pandoc.Div(content,
    pandoc.Attr(div.identifier, {}, { ["custom-style"] = docx_style[kind] }))
  return pandoc.List({ styled })
end

render.html = function(kind, div, content)
  local status = div.attributes["status"]
  div.content = content
  if status then div.classes:insert("status-" .. status) end
  return pandoc.List({ div })
end

------------------------------------------------------------ structure walk

local kind_of, transform

local function walk_blocks(blocks)
  local out = pandoc.List()
  for _, block in ipairs(blocks) do
    if block.t == "Header" and block.level == 1 then
      local mark = nil
      if block.classes:includes("unnumbered") then
        mark = section_mark("")
      elseif block.classes:includes("editor") then
        sec.kind = "editor"
        mark = section_mark(cfg.labels.editor)
      else
        sec.kind = "reviewer"
        reviewer_n = reviewer_n + 1
        comment_n = 0
        mark = section_mark(string.format("%s %d", cfg.labels.reviewer, reviewer_n))
      end
      out:insert(block)
      if mark then out:insert(mark) end
    elseif block.t == "Div" and kind_of(block) then
      out:extend(transform(block))
    elseif block.t == "Div" and block.classes:includes("response-summary") then
      out:insert(block)                       -- placeholder, resolved later
    elseif block.t == "Div" then
      block.content = walk_blocks(block.content)
      out:insert(block)
    else
      out:insert(block)
    end
  end
  return out
end

kind_of = function(div)
  for _, kind in ipairs(KINDS) do
    if div.classes:includes(kind) then return kind end
  end
  return nil
end

transform = function(div)
  local kind = kind_of(div)
  local number = nil

  if kind == "comment" then
    local general = div.classes:includes("general")
    if not general and sec.kind == nil then
      general = true
      if not warned_unsectioned then
        warned_unsectioned = true
        warn("comment outside a reviewer/editor section rendered without a " ..
             "number; start sections with '# Reviewer 1' or '# ... {.editor}'")
      end
    end
    if not general then
      local anchor
      number, anchor = comment_number()
      if div.identifier == "" then div.identifier = anchor end
      last_comment = {
        label   = string.format("%s %s", cfg.labels.comment, number),
        anchor  = div.identifier,
        summary = div.attributes["summary"],
        status  = nil,
      }
      summary_items:insert(last_comment)
    end
  elseif kind == "reply" then
    local status = div.attributes["status"]
    if status and last_comment then last_comment.status = status end
    if status == "todo" then
      todo_items:insert(last_comment and last_comment.label or "(unnumbered comment)")
    end
  end

  local content = walk_blocks(div.content)
  content:insert(1,
    label_block(kind, label_text(kind, div, number), div.attributes["status"]))

  return render[target](kind, div, content)
end

------------------------------------------------------- summary of changes

local function status_suffix(status)
  if status == "done" then
    local inlines = pandoc.Inlines({ pandoc.Space() })
    inlines:extend(text_inlines(cfg.status.done))
    return inlines
  elseif status == "partial" or status == "declined" or status == "todo" then
    return { pandoc.Space(), pandoc.Str("(" .. cfg.status[status] .. ")") }
  end
  return {}
end

local function build_summary()
  local items = pandoc.List()
  for _, item in ipairs(summary_items) do
    if item.summary then
      local inlines = pandoc.Inlines({
        pandoc.Link({ pandoc.Str(item.label) }, "#" .. item.anchor) })
      inlines:extend(status_suffix(item.status))
      inlines:insert(pandoc.Str(" — "))
      inlines:extend(pandoc.utils.blocks_to_inlines(
        pandoc.read(item.summary).blocks))
      items:insert({ pandoc.Plain(inlines) })
    end
  end
  if #items == 0 then
    warn("a summary was requested ('.response-summary' div or " ..
         "'response-letter: summary: true'), but no comment carries a " ..
         "summary=\"...\" attribute; the summary section was dropped")
    return pandoc.List()
  end
  local blocks = pandoc.List()
  blocks:insert(pandoc.Para({ pandoc.Strong({ pandoc.Str(cfg.labels.summary) }) }))
  blocks:insert(pandoc.BulletList(items))
  if target == "docx" then
    return pandoc.List({ pandoc.Div(blocks, pandoc.Attr("", {}, {})) })
  end
  return pandoc.List({ pandoc.Div(blocks, pandoc.Attr("", { "response-summary" })) })
end

local function resolve_summaries(blocks)
  local found = false
  blocks = pandoc.Blocks(blocks):walk({
    Div = function(div)
      if div.classes:includes("response-summary") then
        found = true
        if cfg.letter.summary == false then return {} end
        return build_summary()
      end
    end,
  })
  -- summary: true with no explicit placeholder: insert the summary right
  -- before the first section heading, i.e. after the opening prose.
  if cfg.letter.summary == true and not found then
    local built = build_summary()
    local at = #blocks + 1
    for i, b in ipairs(blocks) do
      if b.t == "Header" and b.level == 1 then
        at = i
        break
      end
    end
    for j = #built, 1, -1 do
      blocks:insert(at, built[j])
    end
  end
  return blocks
end

------------------------------------------------------ letterhead & closing

local function inlines_of(parts)
  local inlines = pandoc.Inlines({})
  for i, part in ipairs(parts) do
    if i > 1 then inlines:insert(pandoc.Str(" ")) end
    if type(part) == "string" then
      inlines:insert(pandoc.Str(part))
    else
      inlines:insert(part)
    end
  end
  return inlines
end

local function letterhead_blocks()
  local L = cfg.letter
  local blocks = pandoc.List()

  -- Re: Manuscript JML-2026-0123 — "A Tutorial on Things"
  if L.id or L.mtitle then
    local parts = { cfg.labels.re }
    if L.id then
      table.insert(parts, cfg.labels.quoted)
      table.insert(parts, L.id)
    end
    local inlines = inlines_of(parts)
    if L.mtitle then
      if L.id then
        inlines:insert(pandoc.Str(" — "))
      else
        inlines:insert(pandoc.Str(" "))
      end
      inlines:insert(pandoc.Quoted("DoubleQuote", { pandoc.Str(L.mtitle) }))
    end
    blocks:insert(pandoc.Para({ pandoc.Strong(inlines) }))
  end

  -- Journal of Memory and Language · Revision 1
  if L.journal or L.round then
    local inlines = pandoc.Inlines({})
    if L.journal then
      inlines:insert(pandoc.Emph({ pandoc.Str(L.journal) }))
    end
    if L.round then
      if L.journal then inlines:insert(pandoc.Str(" · ")) end
      inlines:insert(pandoc.Str(cfg.labels.round .. " " .. L.round))
    end
    blocks:insert(pandoc.Para(inlines))
  end

  -- Dr. Jane Simpson, Action Editor
  if L.editor_name then
    local inlines = pandoc.Inlines({ pandoc.Str(L.editor_name) })
    if L.editor_role then
      inlines:insert(pandoc.Str(", " .. L.editor_role))
    end
    blocks:insert(pandoc.Para(inlines))
  end

  -- Dear Dr. Simpson,
  local salutation = L.salutation
  if not salutation and L.editor_name then
    salutation = string.format("%s %s,", cfg.labels.dear, L.editor_name)
  end
  if salutation then
    blocks:insert(pandoc.Para({ pandoc.Str(salutation) }))
  end

  if #blocks == 0 then return nil end
  return pandoc.Div(blocks, pandoc.Attr("", { "letterhead" }))
end

local function signature_blocks()
  local sig = cfg.signature
  if sig == nil then return pandoc.List() end
  if pandoc.utils.type(sig) == "Blocks" or (sig.t == "MetaBlocks") then
    local blocks = pandoc.Blocks({})
    for _, b in ipairs(sig) do blocks:insert(b) end
    -- Lines of the signature stack as written in the YAML.
    return blocks:walk({
      SoftBreak = function() return pandoc.LineBreak() end,
    })
  end
  -- Inline or string value: stack lines with hard line breaks.
  local inlines = pandoc.Inlines({})
  if type(sig) == "string" then
    local first = true
    for line in sig:gmatch("[^\n]+") do
      if not first then inlines:insert(pandoc.LineBreak()) end
      inlines:insert(pandoc.Str(line))
      first = false
    end
  else
    for _, inline in ipairs(sig) do
      if inline.t == "SoftBreak" then
        inlines:insert(pandoc.LineBreak())
      else
        inlines:insert(inline)
      end
    end
  end
  return pandoc.List({ pandoc.Para(inlines) })
end

local function closing_blocks()
  local L = cfg.letter
  local blocks = pandoc.List()
  if L.closing then
    blocks:insert(pandoc.Para({ pandoc.Str(L.closing) }))
  end
  blocks:extend(signature_blocks())
  if #blocks == 0 then return nil end
  return pandoc.Div(blocks, pandoc.Attr("", { "closing" }))
end

------------------------------------------------------------ theme injection

local function hex6(color)
  local h = color:gsub("#", ""):upper()
  if not h:match("^%x%x%x%x%x%x$") then return nil end
  return h
end

-- Mix `frac` of the color into white (for HTML background tints).
local function tint(hex, frac)
  local out = {}
  for pair in hex:gmatch("%x%x") do
    local c = tonumber(pair, 16)
    table.insert(out, string.format("%02X",
      math.floor(c * frac + 255 * (1 - frac) + 0.5)))
  end
  return table.concat(out)
end

local function inject_theme()
  if not (quarto and quarto.doc) then return end

  if target == "latex" then
    local lines = pandoc.List()
    for _, kind in ipairs(KINDS) do
      local h = hex6(cfg.colors[kind])
      if h then
        lines:insert(string.format("\\definecolor{rl%s}{HTML}{%s}",
          kind:gsub("^%l", string.upper), h))
      else
        warn(string.format("invalid color '%s' for '%s' (expected #RRGGBB)",
          cfg.colors[kind], kind))
      end
    end
    if cfg.running_header then
      local short = cfg.shorttitle and latex_escape(cfg.shorttitle) or ""
      lines:insert("\\usepackage{fancyhdr}")
      lines:insert("\\setlength{\\headheight}{15pt}")
      lines:insert("\\pagestyle{fancy}")
      lines:insert("\\fancyhf{}")
      lines:insert("\\fancyhead[L]{\\small\\itshape " .. short .. "}")
      lines:insert("\\fancyhead[R]{\\small\\rightmark}")
      lines:insert("\\fancyfoot[C]{\\small\\thepage}")
      lines:insert("\\renewcommand{\\headrulewidth}{0.2pt}")
      -- fancyhdr >= 4 auto-marks sections with \markboth{title}{}, which
      -- would clear the right mark this filter sets after each heading;
      -- the filter's \markright calls are the only marks we want.
      lines:insert("\\renewcommand{\\sectionmark}[1]{}")
      lines:insert("\\renewcommand{\\subsectionmark}[1]{}")
    end
    quarto.doc.include_text("in-header", table.concat(lines, "\n"))
  elseif target == "html" then
    local vars = pandoc.List()
    for _, kind in ipairs(KINDS) do
      local h = hex6(cfg.colors[kind])
      if h then
        vars:insert(string.format("  --rl-%s: #%s;", kind, h))
        vars:insert(string.format("  --rl-%s-bg: #%s;", kind,
          tint(h, cfg.tints[kind])))
      end
    end
    quarto.doc.include_text("in-header",
      "<style>\n:root {\n" .. table.concat(vars, "\n") .. "\n}\n</style>")
  elseif target == "typst" then
    -- No static Typst asset ships (unlike header.tex / .css); the colored box
    -- functions are defined here so the palette stays the single source of
    -- truth. Tints use .lighten(100-N%) to match header.tex's rl<Kind>!N mix.
    local lines = pandoc.List()
    for _, kind in ipairs(KINDS) do
      local h = hex6(cfg.colors[kind])
      if h then
        lines:insert(string.format('#let rl-%s-color = rgb("#%s")', kind, h))
      else
        warn(string.format("invalid color '%s' for '%s' (expected #RRGGBB)",
          cfg.colors[kind], kind))
      end
    end
    lines:insert([==[#let rl-comment(body) = block(
  fill: rl-comment-color.lighten(95%), stroke: (left: 2.5pt + rl-comment-color),
  inset: (left: 7pt, right: 7pt, top: 6pt, bottom: 6pt),
  radius: 1pt, above: 10pt, below: 10pt, width: 100%,
)[#text(style: "italic")[#body]]]==])
    lines:insert([==[#let rl-reply(body) = block(
  fill: rl-reply-color.lighten(96%), stroke: (left: 2.5pt + rl-reply-color),
  inset: (left: 7pt, right: 7pt, top: 6pt, bottom: 6pt),
  radius: 1pt, above: 10pt, below: 10pt, width: 100%,
)[#body]]==])
    lines:insert([==[#let rl-changes(body) = block(
  fill: rl-changes-color.lighten(88%), stroke: (left: 2.5pt + rl-changes-color),
  inset: (left: 7pt, right: 7pt, top: 6pt, bottom: 6pt),
  radius: 1pt, above: 10pt, below: 10pt, width: 100%,
)[#body]]==])
    lines:insert([==[#let rl-quoted(body) = block(
  fill: rl-quoted-color.lighten(92%), stroke: (left: 3pt + rl-quoted-color),
  inset: (left: 9pt, right: 7pt, top: 6pt, bottom: 6pt),
  radius: 0pt, above: 8pt, below: 8pt, width: 100%,
)[#body]]==])
    -- Color links/refs/citations to match the LaTeX PDF, which sets colorlinks
    -- with linkcolor/citecolor: MidnightBlue and urlcolor: Teal (Typst renders
    -- them black by default). URLs (string dest) get Teal; internal links,
    -- cross-references, and citations get MidnightBlue.
    lines:insert('#show ref: set text(fill: rgb("#2D2F92"))')
    lines:insert('#show cite: set text(fill: rgb("#2D2F92"))')
    lines:insert([==[#show link: it => text(
  fill: if type(it.dest) == str { rgb("#008080") } else { rgb("#2D2F92") }, it)]==])
    if cfg.running_header then
      lines:insert(string.format('#let rl-short = "%s"',
        typst_str(cfg.shorttitle or "")))
      -- Page header: shorttitle (left) + current section (right), skipped on
      -- the title page. The section text comes from the <rl-mark> markers
      -- section_mark() drops after each heading.
      lines:insert([==[#set page(header: context {
  let n = here().page()
  if n > 1 {
    // last marker on or before this page (matches LaTeX \rightmark / botmark);
    // .before(here()) would miss a section heading at the very top of the page.
    let m = query(<rl-mark>).filter(it => it.location().page() <= n)
    let sec = if m.len() > 0 { m.last().value } else { "" }
    set text(size: 9pt, style: "italic")
    box(width: 100%, stroke: (bottom: 0.2pt), inset: (bottom: 3pt))[
      #rl-short #h(1fr) #sec
    ]
  }
})]==])
    end
    quarto.doc.include_text("in-header", table.concat(lines, "\n"))
  end
end

---------------------------------------------------------------------- main

function Pandoc(doc)
  read_meta(doc.meta)

  local blocks = walk_blocks(doc.blocks)
  blocks = resolve_summaries(blocks)

  local letterhead = letterhead_blocks()
  if letterhead then blocks:insert(1, letterhead) end

  local closing = closing_blocks()
  if closing then blocks:insert(closing) end

  inject_theme()

  if #todo_items > 0 then
    warn(string.format("%d reply/replies still marked TODO: %s",
      #todo_items, table.concat(todo_items, ", ")))
  end

  doc.blocks = blocks
  return doc
end
