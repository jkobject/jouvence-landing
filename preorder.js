const preorderForm = document.querySelector("[data-preorder-form]");
const responseFrame = document.querySelector("[data-response-frame]");
const formStatus = document.querySelector("[data-form-status]");
const submitLabel = document.querySelector("[data-submit-label]");
const submitButton = preorderForm?.querySelector('button[type="submit"]');

let formWasSubmitted = false;

const formMessages = {
  en: {
    success: "Thank you — your non-binding pre-order has been recorded. We will contact you before any payment.",
    recorded: "Pre-order recorded",
    submitting: "Recording your pre-order…",
  },
  fr: {
    success: "Merci — votre précommande sans engagement a bien été enregistrée. Nous vous contacterons avant tout paiement.",
    recorded: "Précommande enregistrée",
    submitting: "Enregistrement de votre précommande…",
  },
  es: {
    success: "Gracias — tu prepedido no vinculante se ha registrado. Te contactaremos antes de cualquier pago.",
    recorded: "Prepedido registrado",
    submitting: "Registrando tu prepedido…",
  },
};

function getFormMessages() {
  return formMessages[document.documentElement.lang] || formMessages.en;
}

function showSuccessfulReservation() {
  if (!formWasSubmitted) return;

  const messages = getFormMessages();
  formStatus.textContent = messages.success;
  submitLabel.textContent = messages.recorded;
  submitButton.disabled = false;
  formWasSubmitted = false;
  preorderForm.reset();
}

function markReservationAsSubmitting() {
  formWasSubmitted = true;
  submitButton.disabled = true;
  submitLabel.textContent = getFormMessages().submitting;
  formStatus.textContent = "";
}

preorderForm?.addEventListener("submit", markReservationAsSubmitting);
responseFrame?.addEventListener("load", showSuccessfulReservation);
