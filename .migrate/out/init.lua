-- ==============================================================================
-- Neovim 配置 (迁移自 Intel Mac 源机，路径已适配 M2 / Homebrew)
-- 注意: 先 source ~/.vimrc，主题 / 透明背景 / 快捷键 / 极简状态栏都在那里定义。
-- Neovim 使用内置 LSP + 内置补全，不加载 Coc / vim-plug。
-- ==============================================================================
vim.cmd.source(vim.fn.expand("~/.vimrc"))

-- 解析可执行文件的实际路径: 优先 Neovim 自身 PATH (vim.fn.exepath)，
-- 其次回退到本机常见安装位置，避免写死源机的 /opt/local/bin (MacPorts)。
local function resolve_exe(name, fallbacks)
	local p = vim.fn.exepath(name)
	if p ~= "" then
		return p
	end
	for _, cand in ipairs(fallbacks or {}) do
		if vim.fn.executable(cand) == 1 then
			return cand
		end
	end
	return name
end

local clangd = resolve_exe("clangd", { "/opt/homebrew/opt/llvm/bin/clangd", "/usr/bin/clangd" })
local python3 = resolve_exe("python3", { "/opt/homebrew/bin/python3", "/usr/bin/python3" })
local pyright_server = resolve_exe("pyright-langserver", { "/opt/homebrew/bin/pyright-langserver" })

vim.lsp.config("clangd", {
	cmd = { clangd },
	filetypes = { "c", "cpp", "objc", "objcpp" },
	root_markers = { ".clangd", "compile_commands.json", "compile_flags.txt", ".git" },
})

vim.lsp.config("pyright", {
	cmd = { pyright_server, "--stdio" },
	filetypes = { "python" },
	root_markers = { "pyrightconfig.json", "pyproject.toml", "setup.cfg", "setup.py", "requirements.txt", "Pipfile", ".git" },
	settings = {
		python = { pythonPath = python3 },
	},
})

-- 只启用本机确实可执行的 LSP，缺失时给出提示而不是静默失败
if vim.fn.executable(clangd) == 1 then
	vim.lsp.enable("clangd")
else
	vim.notify("clangd 未找到，C/C++ LSP 不可用", vim.log.levels.WARN)
end
if vim.fn.executable(pyright_server) == 1 then
	vim.lsp.enable("pyright")
else
	vim.notify("pyright-langserver 未找到，Python LSP 不可用", vim.log.levels.WARN)
end

vim.o.completeopt = "menu,menuone,noselect"
vim.api.nvim_create_autocmd("LspAttach", {
	callback = function(event)
		local client = vim.lsp.get_client_by_id(event.data.client_id)
		if client and client:supports_method("textDocument/completion") then
			vim.lsp.completion.enable(true, client.id, event.buf, { autotrigger = true })
		end
	end,
})

vim.keymap.set("i", "<Tab>", function()
	return vim.fn.pumvisible() == 1 and "<C-n>" or "<Tab>"
end, { expr = true, silent = true })
vim.keymap.set("i", "<S-Tab>", function()
	return vim.fn.pumvisible() == 1 and "<C-p>" or "<C-h>"
end, { expr = true, silent = true })
vim.keymap.set("i", "<CR>", function()
	return vim.fn.pumvisible() == 1 and "<C-y>" or "<CR>"
end, { expr = true, silent = true })
vim.keymap.set("i", "<C-Space>", "<C-x><C-o>", { silent = true })

vim.keymap.set("n", "K", vim.lsp.buf.hover, { desc = "LSP hover" })
vim.keymap.set("n", "KK", vim.lsp.buf.hover, { desc = "LSP hover" })
vim.keymap.set("n", "gd", vim.lsp.buf.definition, { desc = "LSP definition" })
vim.keymap.set("n", "gy", vim.lsp.buf.type_definition, { desc = "LSP type definition" })
vim.keymap.set("n", "gi", vim.lsp.buf.implementation, { desc = "LSP implementation" })
vim.keymap.set("n", "gr", vim.lsp.buf.references, { desc = "LSP references" })
vim.keymap.set("n", "[g", vim.diagnostic.goto_prev, { desc = "Previous diagnostic" })
vim.keymap.set("n", "]g", vim.diagnostic.goto_next, { desc = "Next diagnostic" })
vim.keymap.set("n", "<leader>rn", vim.lsp.buf.rename, { desc = "LSP rename" })
vim.keymap.set("n", "<leader>f", function()
	vim.lsp.buf.format({ async = true })
end, { desc = "LSP format" })
vim.keymap.set("n", "<leader>qf", vim.lsp.buf.code_action, { desc = "LSP code action" })
vim.keymap.set("n", "<leader>ac", vim.lsp.buf.code_action, { desc = "LSP code action" })
vim.keymap.set("n", "<leader>a", function()
	vim.diagnostic.setqflist({ open = true })
end, { desc = "Open diagnostics" })
vim.keymap.set("n", "<leader>o", vim.lsp.buf.document_symbol, { desc = "Document symbols" })

vim.api.nvim_create_user_command("Format", function()
	vim.lsp.buf.format({ async = true })
end, {})
