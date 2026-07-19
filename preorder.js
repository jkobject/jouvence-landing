const preorderForm = document.querySelector("[data-preorder-form]");
const responseFrame = document.querySelector("[data-response-frame]");
const formStatus = document.querySelector("[data-form-status]");
const submitLabel = document.querySelector("[data-submit-label]");
const submitButton = preorderForm?.querySelector('button[type="submit"]');

let formWasSubmitted = false;

function showSuccessfulReservation() {
  if (!formWasSubmitted) return;

  formStatus.textContent = "Thank you — your non-binding pre-order has been recorded. We will contact you before any payment.";
  submitLabel.textContent = "Pre-order recorded";
  submitButton.disabled = false;
  formWasSubmitted = false;
  preorderForm.reset();
}

function markReservationAsSubmitting() {
  formWasSubmitted = true;
  submitButton.disabled = true;
  submitLabel.textContent = "Recording your pre-order…";
  formStatus.textContent = "";
}

preorderForm?.addEventListener("submit", markReservationAsSubmitting);
responseFrame?.addEventListener("load", showSuccessfulReservation);
