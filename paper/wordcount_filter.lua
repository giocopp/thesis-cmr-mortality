function Pandoc(doc)
  local body_words = 0
  local counting = true

  for _, block in ipairs(doc.blocks) do
    if block.t == "Header" then
      local text = pandoc.utils.stringify(block):lower()
      if text:match("^references") or text:match("^bibliography") or text:match("^works cited") then
        counting = false
      end
    elseif block.t == "Div" and block.identifier == "refs" then
      counting = false
    elseif block.t == "RawBlock" and (block.format == "latex" or block.format == "tex") then
      if block.text:match("\\appendix") then
        counting = false
      end
    end

    if counting then
      pandoc.walk_block(block, {
        Str = function(el)
          body_words = body_words + 1
        end
      })
    end
  end

  local count = tostring(body_words)
  for i, block in ipairs(doc.blocks) do
    if block.t == "RawBlock" and (block.format == "latex" or block.format == "tex") then
      block.text = string.gsub(block.text, "@WORDCOUNT_BODY@", count)
    end
  end

  return doc
end
