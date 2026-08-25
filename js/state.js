/* Minimal sessionStorage-backed state shared across the CholeraFlow screens.
   Only metadata is persisted here — raw File objects can't survive JSON
   serialisation, so each screen that needs file bytes keeps them in memory
   for its own lifetime. */
(function () {
  var STORAGE_KEY = 'choleraflow_state';

  function defaultState() {
    return {
      maxStepReached: 1,
      upload: null,
    };
  }

  function load() {
    try {
      var raw = sessionStorage.getItem(STORAGE_KEY);
      if (!raw) return defaultState();
      var parsed = JSON.parse(raw);
      return Object.assign(defaultState(), parsed);
    } catch (e) {
      return defaultState();
    }
  }

  function save(state) {
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify(state));
  }

  window.CholeraFlowState = {
    get: load,
    update: function (patch) {
      var next = Object.assign(load(), patch);
      save(next);
      return next;
    },
    clear: function () {
      sessionStorage.removeItem(STORAGE_KEY);
    },
  };
})();
