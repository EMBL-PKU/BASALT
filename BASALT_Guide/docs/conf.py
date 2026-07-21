"""Sphinx configuration for the BASALT documentation."""

from __future__ import annotations

import os


project = "BASALT"
author = "PKU-EMBL Laboratory"
copyright = "2023–2026, PKU-EMBL Laboratory"

extensions = [
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
]

source_suffix = {".md": "markdown"}
master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "linkify",
    "substitution",
]
myst_heading_anchors = 3

html_theme = "sphinx_book_theme"
html_title = "BASALT"
html_logo = "img/pku-embl-logo.png"
html_favicon = "img/pku-embl-logo.png"
html_static_path = ["_static"]
html_css_files = ["css/basalt.css"]
html_show_sphinx = False
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")

html_theme_options = {
    "repository_url": "https://github.com/PKU-EMBL/BASALT",
    "repository_branch": "master",
    "path_to_docs": "BASALT_Guide/docs",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": True,
    "use_source_button": True,
    "use_download_button": True,
    "home_page_in_toc": True,
    "show_navbar_depth": 1,
    "navigation_with_keys": True,
    "logo": {
        "text": "BASALT",
        "image_light": "img/pku-embl-logo.png",
        "image_dark": "img/pku-embl-logo.png",
    },
}

html_context = {
    "display_github": True,
    "github_user": "PKU-EMBL",
    "github_repo": "BASALT",
    "github_version": "master",
    "conf_py_path": "/BASALT_Guide/docs/",
}

copybutton_prompt_text = r"^(\$ |>>> |\.\.\. )"
copybutton_prompt_is_regexp = True
linkcheck_report_timeouts_as_broken = True
# Figshare resolves this DOI in normal browsers but returns HTTP 403 to
# Sphinx's automated link checker.
linkcheck_ignore = [r"https://doi\.org/10\.6084/m9\.figshare\.22323424"]
