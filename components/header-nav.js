/* Renders the sticky header + 5-step nav shared by every CholeraFlow screen.
   Steps beyond what the user has reached are rendered inert (no page exists
   to send them to yet / they haven't unlocked it). */
(function () {
  var STEP_LABELS = ['Upload & validate', 'Select analysis', 'Configure', 'Run & monitor', 'Results'];
  var STEP_PAGES = ['01-upload.html', 'screen2_select_mode.html', '03-configure.html', '04-run.html', '05-results.html'];

  function renderHeader(container, currentStep) {
    var state = window.CholeraFlowState.get();
    var maxReached = Math.max(state.maxStepReached || 1, currentStep);

    var header = document.createElement('header');
    header.className = 'cf-header';

    var brand = document.createElement('div');
    brand.className = 'cf-header__brand';
    brand.innerHTML =
      '<span class="cf-logo">C</span>' +
      '<span class="cf-header__title">CholeraFlow</span>' +
      '<span class="cf-pill cf-pill--tag cf-pill--teal">Public health build</span>';

    var nav = document.createElement('nav');
    nav.className = 'cf-nav';
    STEP_LABELS.forEach(function (label, i) {
      var n = i + 1;
      var isActive = n === currentStep;
      var isDone = n < currentStep;
      var reachable = n <= maxReached;

      var btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'cf-nav__step';
      btn.title = label;
      if (!reachable || isActive) btn.disabled = true;
      if (reachable && !isActive) {
        btn.addEventListener('click', function () { window.location.href = STEP_PAGES[i]; });
      }

      var circle = document.createElement('span');
      circle.className = 'cf-nav__circle' + (isActive ? ' cf-nav__circle--active' : isDone ? ' cf-nav__circle--done' : '');
      circle.textContent = isDone ? '✓' : String(n);

      var labelEl = document.createElement('span');
      labelEl.className = 'cf-nav__label' + (isActive ? ' cf-nav__label--active' : '');
      labelEl.textContent = label;

      btn.appendChild(circle);
      btn.appendChild(labelEl);
      nav.appendChild(btn);
    });

    var right = document.createElement('div');
    right.className = 'cf-header__right';
    right.innerHTML =
      '<span class="cf-header__help" title="Open the help centre">Help</span>' +
      '<span class="cf-avatar">AL</span>';

    header.appendChild(brand);
    header.appendChild(nav);
    header.appendChild(right);

    container.replaceChildren(header);
  }

  window.CholeraFlowHeader = { render: renderHeader };
})();
