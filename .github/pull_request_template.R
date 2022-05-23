**What changes are proposed in this pull request?**


**If there is an GitHub issue associated with this pull request, please provide link.**


--------------------------------------------------------------------------------

  Reviewer Checklist (if item does not apply, mark is as complete)

- [ ] PR branch has pulled the most recent updates from main branch. Ensure the pull request branch and your local version match and both have the latest updates from the main branch.
- [ ] If a new function was added, function included in `_pkgdown.yml`
- [ ] If a bug was fixed, a unit test was added for the bug check
- [ ] Run `pkgdown::build_site()`. Check the R console for errors, and review the rendered website.
- [ ] Code coverage is suitable for any new functions/features. Review coverage with `withr::with_envvar(new = c("NOT_CRAN" = "true"), covr::report())`. Begin in a fresh R session without any packages loaded.
- [ ] R CMD Check runs without errors, warnings, and notes
- [ ] `usethis::use_spell_check()` runs with no spelling errors in documentation

When the branch is ready to be merged into master:
 - [ ] Update `NEWS.md` with the changes from this pull request under the heading "`# cbioportalR (development version)`". If there is an issue associated with the pull request, reference it in parentheses at the end update (see `NEWS.md` for examples).
- [ ] Run `codemetar::write_codemeta()`
- [ ] Run `usethis::use_spell_check()` again
- [ ] Approve Pull Request
- [ ] Merge the PR
