# EasyMath Design System Skill

## Overview
This skill provides the complete design system for EasyMath interactive demos. Use this when creating new demo pages or components to ensure consistency with the established design language.

## Quick Start
When creating a new demo page:
1. Copy the theme structure from `demo/src/styles/themes.css`
2. Import global styles from `demo/src/styles/global.css`
3. Use the layout template and component patterns below
4. Follow the spacing, typography, and color guidelines

## Design Tokens

### Colors
All colors use CSS variables for theme support. Key variables:
- `--bg-primary`, `--bg-secondary`, `--bg-tertiary`, `--bg-quaternary`
- `--text-primary`, `--text-secondary`, `--text-muted`
- `--accent-primary`, `--accent-secondary`, `--accent-primary-alpha`
- `--border-color`, `--border-subtle`

### Spacing
- `--space-xs`: 4px
- `--space-sm`: 8px
- `--space-md`: 12px
- `--space-lg`: 16px
- `--space-xl`: 24px
- `--space-2xl`: 32px

### Typography
- `--text-xs`: 0.75rem
- `--text-sm`: 0.875rem
- `--text-base`: 1rem
- `--text-md`: 1.125rem
- `--text-lg`: 1.25rem
- `--text-xl`: 1.5rem

### Border Radius
- `--radius-sm`: 4px
- `--radius-md`: 6px
- `--radius-lg`: 12px
- `--radius-full`: 9999px

### Shadows
- `--shadow-sm`: 0 1px 2px rgba(0,0,0,0.1)
- `--shadow-md`: 0 2px 8px rgba(0,0,0,0.2)
- `--shadow-lg`: 0 4px 16px rgba(0,0,0,0.3)
- `--shadow-xl`: 0 8px 32px rgba(0,0,0,0.4)
- `--shadow-glow`: 0 0 12px var(--accent-primary-alpha)
- `--shadow-glow-strong`: 0 0 20px var(--accent-primary-alpha)

### Transitions
- `--transition-fast`: 150ms ease
- `--transition-base`: 200ms ease
- `--transition-slow`: 300ms ease

### Z-Index Layers
- `--z-base`: 1
- `--z-resizer`: 10
- `--z-header`: 50
- `--z-dropdown`: 200
- `--z-modal`: 1000

## Component Patterns

### 1. Button Styles

**Primary Button**:
```css
.btn-primary {
  padding: var(--space-sm) var(--space-lg);
  background: var(--accent-primary-alpha);
  border: 1px solid var(--accent-primary);
  border-radius: var(--radius-md);
  color: var(--accent-primary);
  cursor: pointer;
  transition: all var(--transition-fast);
}

.btn-primary:hover {
  background: rgba(var(--accent-primary), 0.25);
  transform: translateY(-1px);
  box-shadow: var(--shadow-glow);
}
```

**Secondary Button**:
```css
.btn {
  padding: var(--space-sm) var(--space-lg);
  background: var(--bg-tertiary);
  border: 1px solid var(--border-color);
  border-radius: var(--radius-md);
  color: var(--text-secondary);
  cursor: pointer;
  transition: all var(--transition-fast);
}

.btn:hover {
  background: var(--bg-quaternary);
  color: var(--text-primary);
}
```

### 2. Slider Control
```html
<div class="slider-group">
  <div class="row">
    <div class="row-header">
      <label>Parameter Name</label>
      <span id="param-v">0.0</span>
    </div>
    <input type="range" id="param" min="-3" max="3" step="0.1" value="0">
  </div>
</div>
```

```css
.row {
  display: flex;
  flex-direction: column;
  gap: var(--space-sm);
  margin-bottom: var(--space-lg);
}

.row-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
}

input[type=range] {
  flex: 1;
  accent-color: var(--accent-primary);
  cursor: pointer;
}
```

### 3. Collapsible Panel
```html
<details class="cs" open>
  <summary class="cs-summary">
    Panel Title
    <svg class="cs-chevron">...</svg>
  </summary>
  <div class="cs-content">
    Panel content here
  </div>
</details>
```

```css
.cs {
  margin-bottom: var(--space-lg);
  border: 1px solid var(--border-subtle);
  border-radius: var(--radius-md);
  background: rgba(255, 255, 255, 0.01);
}

.cs-summary {
  padding: var(--space-sm) var(--space-md);
  cursor: pointer;
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.cs-content {
  padding: var(--space-md);
  animation: slideDown 0.2s ease;
}

@keyframes slideDown {
  from { opacity: 0; transform: translateY(-4px); }
  to { opacity: 1; transform: translateY(0); }
}
```

### 4. Tab Switcher
```html
<div id="tab-switcher">
  <button class="tab-btn active" data-tab="tab1">
    <span>Tab 1</span>
  </button>
  <button class="tab-btn" data-tab="tab2">
    <span>Tab 2</span>
  </button>
</div>
```

```css
#tab-switcher {
  display: inline-flex;
  background: var(--bg-secondary);
  border: 1px solid var(--border-color);
  border-radius: var(--radius-md);
  padding: 3px;
  gap: 2px;
}

.tab-btn {
  padding: var(--space-sm) 18px;
  border: none;
  background: transparent;
  color: var(--text-secondary);
  border-radius: var(--radius-sm);
  transition: all var(--transition-fast);
}

.tab-btn.active {
  background: var(--bg-tertiary);
  color: var(--accent-primary);
  box-shadow: var(--shadow-glow);
}
```

### 5. Icon Button with Badge
```html
<button class="icon-btn">
  <svg>...</svg>
  <span class="badge">Default</span>
</button>
```

### 6. Dropdown Menu
```html
<div class="dropdown">
  <button class="dropdown-btn" aria-haspopup="true" aria-expanded="false">
    Menu
    <svg class="chevron">...</svg>
  </button>
  <div class="dropdown-menu" role="menu" aria-hidden="true">
    <button class="dropdown-item" role="menuitem">Option 1</button>
    <button class="dropdown-item" role="menuitem">Option 2</button>
  </div>
</div>
```

## Layout Patterns

### Three-Panel Layout (Desktop)
```css
.app {
  display: grid;
  grid-template-columns: var(--sidebar-width, 280px) 4px 1fr;
  grid-template-rows: var(--header-height) 1fr auto var(--bottom-height);
  height: 100vh;
  overflow: hidden;
}

.header { grid-column: 1 / -1; grid-row: 1; }
.sidebar { grid-column: 1; grid-row: 2 / -1; }
.canvas { grid-column: 3; grid-row: 2; }
.panel-bottom { grid-column: 3; grid-row: 4; }
```

### Responsive Breakpoints
```css
/* Desktop > 1024px */
.app {
  grid-template-columns: var(--sidebar-width, 280px) 4px 1fr;
}

/* Tablet 768px - 1024px */
@media (max-width: 1024px) {
  .app {
    grid-template-columns: 240px 4px 1fr;
  }
}

/* Mobile < 768px */
@media (max-width: 768px) {
  .app {
    grid-template-columns: 1fr;
    grid-template-rows: var(--header-height) minmax(300px, 1fr) auto 200px;
  }
  .sidebar {
    position: fixed;
    transform: translateX(-100%);
  }
}
```

## Theme Implementation

### 1. Define Theme Variables
```css
:root[data-theme="my-theme"] {
  --bg-primary: #0a0a0f;
  --bg-secondary: #111827;
  --accent-primary: #8b7aff;
  --accent-gradient: linear-gradient(135deg, #8b7aff 0%, #5eb3ff 100%);
  /* ... all other variables */
}
```

### 2. Theme Switcher Component
```astro
<ThemeSwitcher />
```

### 3. Persist Theme Preference
```javascript
const savedTheme = localStorage.getItem('theme-preference') || 'purple-blue';
document.documentElement.setAttribute('data-theme', savedTheme);
```

## Common Interactions

### Parameter Feedback Pattern
```javascript
function onParameterChange(paramId, value) {
  // Update display
  document.getElementById(paramId + '-v').textContent = value;

  // Highlight affected components
  highlightAffectedComponents(paramId);

  // Auto-clear highlights after delay
  setTimeout(() => clearHighlights(), 2000);
}
```

### Panel State Persistence
```javascript
// Save on change
details.addEventListener('toggle', () => {
  localStorage.setItem('panel-state', details.open);
});

// Restore on load
const savedState = localStorage.getItem('panel-state');
details.open = savedState === 'true';
```

### Resize Handling
```javascript
let isResizing = false;
let startPos = 0;

resizer.addEventListener('mousedown', (e) => {
  isResizing = true;
  startPos = e.clientX;
  document.body.style.cursor = 'col-resize';
});

document.addEventListener('mousemove', (e) => {
  if (!isResizing) return;
  const delta = e.clientX - startPos;
  updateSize(delta);
});

document.addEventListener('mouseup', () => {
  if (isResizing) {
    isResizing = false;
    document.body.style.cursor = '';
    saveSize();
  }
});
```

## Best Practices

### 1. Accessibility
- Always include ARIA labels on interactive elements
- Support keyboard navigation (Tab, Enter, Escape)
- Provide focus indicators
- Ensure color contrast meets WCAG AA
- Support reduced motion preferences

### 2. Performance
- Use GPU-accelerated animations (transform, opacity)
- Throttle expensive calculations (60fps max)
- Debounce localStorage writes
- Use CSS variables for theme switching (instant)

### 3. Responsive Design
- Design mobile-first
- Use relative units (rem, em, %)
- Test on actual devices
- Ensure touch targets are ≥44px
- Simplify layouts on small screens

### 4. Code Organization
- Keep styles co-located with components
- Use consistent naming conventions
- Document complex interactions
- Separate concerns (structure, style, behavior)

### 5. User Experience
- Provide immediate feedback
- Show clear cause-and-effect relationships
- Remember user preferences
- Use animations to guide attention
- Keep interfaces spacious and breathable

## Resources
- Full design specs: `memory/matrix-demo-design.md`
- Theme definitions: `demo/src/styles/themes.css`
- Component examples: `demo/src/pages/matrix/index.astro`
- Global styles: `demo/src/styles/global.css`
