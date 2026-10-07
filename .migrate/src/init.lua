vim.cmd.source(vim.fn.expand("~/.vimrc"))

local node = vim.fn.expand("~/.local/bin/node")
local pyright_server = vim.fn.expand("~/.config/coc/extensions/node_modules/pyright/langserver.index.js")

vim.lsp.config("clangd", {
	cmd = { "/opt/local/bin/clangd" },
	filetypes = { "c", "cpp", "objc", "objcpp" },
	root_markers = { ".clangd", "compile_commands.json", "compile_flags.txt", ".git" },
})

vim.lsp.config("pyright", {
	cmd = { node, pyright_server, "--stdio" },
	filetypes = { "python" },
	root_markers = { "pyrightconfig.json", "pyproject.toml", "setup.cfg", "setup.py", "requirements.txt", "Pipfile", ".git" },
	settings = {
		python = { pythonPath = "/opt/local/bin/python3" },
	},
})

vim.lsp.enable({ "clangd", "pyright" })

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