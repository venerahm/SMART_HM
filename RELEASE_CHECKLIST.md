# Release Checklist

Use this checklist when creating the publication-ready GitHub repository.

## 1. Validate This Working Repository

From the current working repository:

```bash
Rscript scripts/validate_publication_manifest.R
```

The validation should report:

```text
Publication manifest validation passed.
```

## 2. Create The Clean Repository Snapshot

Choose a destination outside this working repository, then run:

```bash
Rscript scripts/create_publication_snapshot.R /path/to/new/repo \
  --overwrite \
  --init-git \
  --remote-url=https://github.com/<user>/<repo>.git \
  --commit-message="Initial publication-ready repository"
```

The snapshot script copies only files marked for publication in
`PUBLICATION_MANIFEST.tsv`, initializes git on branch `main`, stages the files,
configures `origin` when a remote URL is supplied, and creates the initial
commit when a commit message is supplied. By default, the script refuses to
write the clean repository inside the source working repository, which avoids
accidental nested git repositories.

## 3. Validate The Clean Repository

From the clean repository:

```bash
Rscript scripts/validate_publication_manifest.R
git status --short
git log --oneline --decorate -1
```

Expected state:

- validation passes;
- `git status --short` prints nothing;
- the latest commit is `Initial publication-ready repository`;
- branch is `main`;
- repository size is small enough for ordinary GitHub use.

## 4. Push To GitHub

From the clean repository:

```bash
git push -u origin main
```

If the remote URL was not supplied during snapshot creation, add it first:

```bash
git remote add origin https://github.com/<user>/<repo>.git
git push -u origin main
```

## 5. Public Release Review

Before sharing the repository URL, confirm:

- `README.md` describes the current manuscript workflow;
- `REQUIREMENTS.md` explains R packages and external data inputs;
- selected summary CSV files and figure sources are present;
- raw iteration outputs, local workspaces, logs, and temporary files are absent;
- participant-level input data are absent.
