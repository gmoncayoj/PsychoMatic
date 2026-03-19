system("git --version")
usethis::git_sitrep()

usethis::use_git_config(user.name = "gmoncayoj",
                        user.email = "gamarramoncayoj@gmail.com")

usethis::use_git()
usethis::use_github()

usethis::create_github_token()
gitcreds::gitcreds_set()
