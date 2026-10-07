-- 迁移验证脚本 (由 nvim --headless 执行)
local function out(k, v)
	print(k .. "=" .. tostring(v))
end

out("NVIM_VERSION", vim.fn.hostname() .. "/" .. tostring(vim.version()))
out("COLORSCHEME", tostring(vim.g.colors_name))
out("TERMGUICOLORS", tostring(vim.o.termguicolors))
out("LASTSTATUS", tostring(vim.o.laststatus))
out("CMDHEIGHT", tostring(vim.o.cmdheight))
out("STATUSLINE", tostring(vim.o.statusline))
out("MAPLEADER", tostring(vim.g.mapleader))

-- Coc 必须不存在
out("EXISTS_CocList", vim.fn.exists(":CocList"))
out("EXISTS_CocCommand", vim.fn.exists(":CocCommand"))
out("RT_HAS_COC", tostring(vim.o.runtimepath:find("coc") ~= nil))
out("HAS_COC_LUA", tostring(pcall(require, "coc")))
local coc_cmds = vim.tbl_filter(function(c)
	return c:lower():find("^coc") ~= nil
end, vim.fn.getcompletion("", "command"))
out("COC_COMMANDS", table.concat(coc_cmds, ","))

-- 依赖路径解析
out("EXEPATH_clangd", vim.fn.exepath("clangd"))
out("EXEPATH_python3", vim.fn.exepath("python3"))
out("EXEPATH_pyright-langserver", vim.fn.exepath("pyright-langserver"))
out("EXEC_clangd", vim.fn.executable(vim.fn.exepath("clangd")))
out("EXEC_pyright", vim.fn.executable(vim.fn.exepath("pyright-langserver")))

-- 原生 LSP 实际启动测试
local function probe(ft, file)
	vim.cmd("silent! edit " .. file)
	vim.bo.filetype = ft
	vim.wait(8000, function()
		return #vim.lsp.get_clients({ bufnr = 0 }) > 0
	end, 100)
	local names = {}
	for _, c in ipairs(vim.lsp.get_clients({ bufnr = 0 })) do
		names[#names + 1] = c.name .. "(" .. tostring(c.config.cmd[1]) .. ")"
	end
	out("LSP_" .. ft, #names > 0 and table.concat(names, ",") or "NONE")
	vim.cmd("silent! bwipeout!")
end

probe("c", "/Users/junxiongnie/code/algo/ParseLQCData/.migrate/probe/test.c")
probe("python", "/Users/junxiongnie/code/algo/ParseLQCData/.migrate/probe/test.py")

out("VERIFY_DONE", "yes")
