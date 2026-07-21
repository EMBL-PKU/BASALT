document.addEventListener("DOMContentLoaded", function () {
  const searchInput = document.querySelector(".md-search__input");
  if (!searchInput) return;

  searchInput.addEventListener("focus", function () {
    document.dispatchEvent(new CustomEvent("readthedocs-search-show"));
  });
});
