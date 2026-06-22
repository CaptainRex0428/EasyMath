/* ============================================
   storage.js - Safe localStorage wrapper
   ============================================ */

const PREFIX = 'easymath-';

const safeStorage = {
  /**
   * Get a value from localStorage. Returns the default if missing or unavailable.
   */
  get(key, defaultValue = null) {
    try {
      const raw = localStorage.getItem(PREFIX + key);
      if (raw === null) return defaultValue;
      return JSON.parse(raw);
    } catch (e) {
      // localStorage unavailable (private mode) or value not JSON
      return defaultValue;
    }
  },

  /**
   * Set a value in localStorage. Silently fails if unavailable.
   */
  set(key, value) {
    try {
      localStorage.setItem(PREFIX + key, JSON.stringify(value));
      return true;
    } catch (e) {
      return false;
    }
  },

  /**
   * Remove a key from localStorage.
   */
  remove(key) {
    try {
      localStorage.removeItem(PREFIX + key);
      return true;
    } catch (e) {
      return false;
    }
  }
};

// Export to global scope for use in inline scripts
window.EasyMathStorage = safeStorage;