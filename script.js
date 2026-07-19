const header = document.querySelector("[data-header]");
const menuButton = document.querySelector("[data-menu-button]");
const mobileNavigation = document.querySelector("[data-mobile-nav]");
const revealElements = document.querySelectorAll(".reveal");
const permissionToggles = document.querySelectorAll(".toggle");

function updateHeader() {
  header?.classList.toggle("is-scrolled", window.scrollY > 40);
}

function closeMobileMenu() {
  if (!menuButton || !mobileNavigation || !header) return;
  menuButton.setAttribute("aria-expanded", "false");
  mobileNavigation.setAttribute("aria-hidden", "true");
  mobileNavigation.inert = true;
  mobileNavigation.classList.remove("is-open");
  header.classList.remove("menu-open");
  document.body.style.overflow = "";
}

function toggleMobileMenu() {
  if (!menuButton || !mobileNavigation || !header) return;
  const isOpen = menuButton.getAttribute("aria-expanded") === "true";
  menuButton.setAttribute("aria-expanded", String(!isOpen));
  mobileNavigation.setAttribute("aria-hidden", String(isOpen));
  mobileNavigation.inert = isOpen;
  mobileNavigation.classList.toggle("is-open", !isOpen);
  header.classList.toggle("menu-open", !isOpen);
  document.body.style.overflow = isOpen ? "" : "hidden";
}

function togglePermission(event) {
  const toggle = event.currentTarget;
  const isEnabled = toggle.getAttribute("aria-checked") === "true";
  toggle.setAttribute("aria-checked", String(!isEnabled));
  toggle.classList.toggle("is-on", !isEnabled);
}

function revealHashTarget() {
  if (!window.location.hash) return;

  const target = document.getElementById(window.location.hash.slice(1));
  if (!target) return;

  target.classList.add("is-visible");
  target.querySelectorAll(".reveal").forEach((element) => element.classList.add("is-visible"));
}

const revealObserver = new IntersectionObserver(
  (entries) => {
    entries.forEach((entry) => {
      if (entry.isIntersecting) {
        entry.target.classList.add("is-visible");
        revealObserver.unobserve(entry.target);
      }
    });
  },
  { threshold: 0.12 },
);

revealElements.forEach((element) => revealObserver.observe(element));
revealHashTarget();
permissionToggles.forEach((toggle) => toggle.addEventListener("click", togglePermission));
menuButton?.addEventListener("click", toggleMobileMenu);
mobileNavigation?.querySelectorAll("a").forEach((link) => link.addEventListener("click", closeMobileMenu));
document.addEventListener("keydown", (event) => {
  if (event.key === "Escape" && menuButton?.getAttribute("aria-expanded") === "true") {
    closeMobileMenu();
    menuButton.focus();
  }
});
window.addEventListener("resize", () => {
  if (window.innerWidth > 1050 && menuButton?.getAttribute("aria-expanded") === "true") closeMobileMenu();
});
window.addEventListener("scroll", updateHeader, { passive: true });
window.addEventListener("hashchange", revealHashTarget);
document.querySelectorAll("[data-year]").forEach((element) => {
  element.textContent = new Date().getFullYear();
});
updateHeader();
