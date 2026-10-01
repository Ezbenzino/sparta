# 代码发布清单（git → GitHub → Zenodo）

> **代码可公开获取是 Bioinformatics 与 BIB 的硬性投稿条件**，不是加分项。
> 这一步和数据无关，现在就能做完。

---

## 0. 先清掉我留下的东西（一条命令）

我试着从 Cowork 的桥接侧初始化 git，失败了，原因值得记一笔：

> git 每次写索引都要**先创建 `.git/index.lock`、写完再删掉它**。
> 而挂载到桥接侧的目录**禁止删除文件**（`rm` 一律 Operation not permitted）。
> 于是每条 git 命令都成功地创建了锁、却删不掉，给下一条命令留了个绊子。
> **结论：这个仓库的 git 操作必须在 Windows 本地做，不能走桥接。**

`.git/` 目前是空仓库（**没有任何 commit**），里面留了两个改名后的锁和一个临时文件。
最干净的处理是整个删掉重来：

```powershell
cd D:\sparta
Remove-Item -Recurse -Force .git
```

另外 `downloads\_sparta_src_for_check.tgz` 是我做校验用的打包，可以删（已被 gitignore 覆盖，不删也不会进仓库）。

---

## 1. 填掉三处占位符（提交前必须做）

| 文件 | 要填什么 |
|---|---|
| `LICENSE` | `<YOUR NAME OR INSTITUTION>` — 版权行 |
| `CITATION.cff` | `<SURNAME>` / `<GIVEN NAME>`，可选 ORCID 与单位；`repository-code` 等建库后回填 |
| `README.md` | 末尾无占位符；`docs/manuscript_abstract.md` 的 Availability 段有 `<GitHub URL>` 与 `<Zenodo DOI>`，拿到后回填 |

---

## 2. 建库并首次提交

```powershell
cd D:\sparta
git init
git add -A
git status --short          # 核对下面那张表
git commit -m "SPARTA v2.0: dual spatial barrier framework"
git branch -M main
```

### 提交前的核对表（我已在隔离环境验证过）

我把源码树复制到一个可删除的环境里跑了一遍 `git add -A`，结果应当是：

| 项 | 期望值 |
|---|---|
| 文件数 | **95** |
| 体积合计 | **0.67 MB** |
| 最大的文件 | `docs/manuscript_results.md`（33 KB） |
| `data/` 下 | 只有 `ledger.csv` 与 `admission_audit.csv`（都是元数据，应当入库） |
| `results/` 下 | 只有 `.gitkeep` |
| 原始数据 / 中间产物 / 图件 | **一个都没有** |

**如果 `git status --short` 显示的文件数量级不对（比如上千个、或几百 MB），
先停下来**，多半是 `.gitignore` 没生效或者在错误的目录里 init 了。

我在 `.gitignore` 里补了两条：`downloads/`（原来没盖住，里面有个 38 MB 的
`GSE250636_RAW.tar`，`git add .` 会把它提交进历史，之后再删也洗不掉）和 `*.bak_*`。

---

## 3. 推到 GitHub（公开）

```powershell
gh repo create sparta --public --source=. --remote=origin --push
# 没装 gh 的话：在网页上建空仓库，然后
#   git remote add origin https://github.com/<USER>/sparta.git
#   git push -u origin main
```

推完把 URL 填回 `CITATION.cff` 的 `repository-code`。

---

## 4. 打 tag 并归档到 Zenodo

```powershell
git tag -a v2.0.0 -m "SPARTA v2.0.0"
git push origin v2.0.0
```

然后：

1. 登录 https://zenodo.org，用 GitHub 账号授权
2. 在 **Settings → GitHub** 里把 `sparta` 仓库的开关打开
3. 回到 GitHub，基于 `v2.0.0` 建一个 **Release**（Zenodo 只在建 Release 时归档，
   光打 tag 不触发）
4. Zenodo 会自动生成 DOI，把它填回 `CITATION.cff` 与摘要的 Availability 段

**顺序很重要**：必须先在 Zenodo 里打开开关，再建 Release。反了的话第一个 Release
不会被归档，得再发一版。

---

## 5. 投稿前最后核对

- [ ] `LICENSE` 与 `CITATION.cff` 里没有 `<>` 占位符
- [ ] 仓库是 **public**（很多人建成 private 然后忘了改）
- [ ] README 的 Data availability 段列全了 GEO accession
- [ ] Zenodo DOI 已填进摘要的 Availability 段
- [ ] `git log` 里没有任何数据文件（`git log --stat | Select-String "\.h5ad|\.tar|\.npz"` 应当为空）

最后一条尤其重要：数据一旦进了 git 历史，**删掉文件不会让它从历史里消失**，
只能 rewrite history。投稿前查一次，比投完被编辑问起来强。
