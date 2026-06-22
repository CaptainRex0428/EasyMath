/* ============================================
   resizable-sidebar.js
   Sidebar drag-to-resize with localStorage persistence
   ============================================ */

(function () {
  const STORAGE_KEY = 'sidebar-width';
  const MIN_WIDTH = 200;
  const MAX_WIDTH = 600;
  const DEFAULT_WIDTH = 280;

  function initSidebar() {
    const sidebar = document.getElementById('controls');
    const resizer = document.getElementById('sidebar-resizer');
    if (!sidebar || !resizer) return;

    // Restore saved width
    const saved = window.EasyMathStorage.get(STORAGE_KEY, DEFAULT_WIDTH);
    const validWidth = Math.min(Math.max(saved, MIN_WIDTH), MAX_WIDTH);
    sidebar.style.width = validWidth + 'px';
    document.documentElement.style.setProperty('--sidebar-width', validWidth + 'px');

    let isResizing = false;
    let startX = 0;
    let startWidth = validWidth;

    resizer.addEventListener('mousedown', (e) => {
      isResizing = true;
      startX = e.clientX;
      startWidth = sidebar.getBoundingClientRect().width;
      document.body.style.cursor = 'col-resize';
      document.body.style.userSelect = 'none';
      resizer.classList.add('resizing');
      e.preventDefault();
    });

    document.addEventListener('mousemove', (e) => {
      if (!isResizing) return;
      const delta = e.clientX - startX;
      const newWidth = Math.min(Math.max(startWidth + delta, MIN_WIDTH), MAX_WIDTH);
      sidebar.style.width = newWidth + 'px';
      document.documentElement.style.setProperty('--sidebar-width', newWidth + 'px');

      // Notify resize listeners (for canvas resize)
      window.dispatchEvent(new CustomEvent('easymath:resize'));
    });

    document.addEventListener('mouseup', () => {
      if (!isResizing) return;
      isResizing = false;
      document.body.style.cursor = '';
      document.body.style.userSelect = '';
      resizer.classList.remove('resizing');
      const finalWidth = sidebar.getBoundingClientRect().width;
      window.EasyMathStorage.set(STORAGE_KEY, Math.round(finalWidth));
    });

    // Touch support for tablets
    resizer.addEventListener('touchstart', (e) => {
      if (e.touches.length !== 1) return;
      isResizing = true;
      startX = e.touches[0].clientX;
      startWidth = sidebar.getBoundingClientRect().width;
      resizer.classList.add('resizing');
    }, { passive: true });

    document.addEventListener('touchmove', (e) => {
      if (!isResizing || e.touches.length !== 1) return;
      const delta = e.touches[0].clientX - startX;
      const newWidth = Math.min(Math.max(startWidth + delta, MIN_WIDTH), MAX_WIDTH);
      sidebar.style.width = newWidth + 'px';
      document.documentElement.style.setProperty('--sidebar-width', newWidth + 'px');
      window.dispatchEvent(new CustomEvent('easymath:resize'));
    }, { passive: true });

    document.addEventListener('touchend', () => {
      if (!isResizing) return;
      isResizing = false;
      resizer.classList.remove('resizing');
      const finalWidth = sidebar.getBoundingClientRect().width;
      window.EasyMathStorage.set(STORAGE_KEY, Math.round(finalWidth));
    });
  }

  // Initialize when DOM is ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initSidebar);
  } else {
    initSidebar();
  }
})();