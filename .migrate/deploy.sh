#!/bin/bash
# ==============================================================================
# Neovim / Vim 配置迁移部署脚本 (从 Intel Mac 源机 -> 本机 M2)
# 幂等、可重复执行；覆盖前一律备份到 ~/.config-migrate-backup/<时间戳>/
# ==============================================================================
set -uo pipefail

SRC="$(cd "$(dirname "$0")" && pwd)/out"
STAMP="$(date +%Y%m%d-%H%M%S)"
BK="$HOME/.config-migrate-backup/$STAMP"

echo "==> 源目录: $SRC"
echo "==> 备份目录: $BK"

# ---------- 1. 只安装缺失的依赖 ----------
echo "==> 检查依赖"
missing=()
for c in nvim clangd python3 node pyright-langserver; do
  if command -v "$c" >/dev/null 2>&1; then
    printf '    [OK]      %-20s %s\n' "$c" "$(command -v "$c")"
  else
    printf '    [MISSING] %-20s\n' "$c"
    missing+=("$c")
  fi
done

if ! command -v pyright-langserver >/dev/null 2>&1; then
  echo "==> 安装 pyright (Homebrew)"
  brew install pyright || { echo "!! brew install pyright 失败"; exit 1; }
fi
if ! command -v nvim >/dev/null 2>&1; then
  echo "==> 安装 neovim (Homebrew)"
  brew install neovim || { echo "!! brew install neovim 失败"; exit 1; }
fi

# ---------- 2. 备份现有目标文件 ----------
echo "==> 备份现有配置"
mkdir -p "$BK"
for f in "$HOME/.vimrc" "$HOME/.vim/coc-settings.json" "$HOME/.zshrc"; do
  if [ -e "$f" ]; then
    rel="${f#$HOME/}"
    mkdir -p "$BK/$(dirname "$rel")"
    cp -p "$f" "$BK/$rel"
    echo "    backed up: $rel"
  else
    echo "    (不存在，跳过): ${f#$HOME/}"
  fi
done

# ---------- 3. 部署合并后的配置 ----------
echo "==> 部署配置文件"
mkdir -p "$HOME/.config/nvim" "$HOME/.vim"
install -m 644 "$SRC/init.lua"         "$HOME/.config/nvim/init.lua"
install -m 644 "$SRC/vimrc"            "$HOME/.vimrc"
install -m 644 "$SRC/coc-settings.json" "$HOME/.vim/coc-settings.json"
echo "    ~/.config/nvim/init.lua"
echo "    ~/.vimrc"
echo "    ~/.vim/coc-settings.json"

# ---------- 4. 合并 ~/.zshrc (不整份覆盖，幂等) ----------
echo "==> 合并 ~/.zshrc"
python3 - <<'PYEOF'
import pathlib
p = pathlib.Path.home() / ".zshrc"
if not p.exists():
    print("    ~/.zshrc 不存在，跳过")
    raise SystemExit(0)
s = p.read_text()
marker = "# --- nvim editor migration ---"
endmarker = "# --- end nvim editor migration ---"
if marker in s:
    print("    migration 块已存在，跳过（避免重复）")
    raise SystemExit(0)
block = (
    marker + "\n"
    'export EDITOR="nvim"\n'
    'export VISUAL="nvim"\n'
    'export FRESH_EDITOR="fresh"   # 原 EDITOR 值，保留以便切回 fresh\n'
    'alias vim="nvim"\n'
    'alias vi="nvim"\n'
    + endmarker
)
old = 'export EDITOR="fresh"'
if old in s:
    s = s.replace(old, block, 1)
    print("    已就地替换原有的 export EDITOR=\"fresh\" 行")
else:
    s = s.rstrip("\n") + "\n\n" + block + "\n"
    print("    未找到原行，已追加 migration 块")
p.write_text(s)
PYEOF

echo "==> 完成"
echo "$BK" > "$HOME/.config-migrate-backup/LAST_BACKUP"
