# MatrixDemo Design Specifications

## Design Goals
MatrixDemo is an interactive teaching tool that helps students understand 3D transformations and the graphics rendering pipeline. The design must clearly show the relationship between:
1. User-adjusted parameters (sliders)
2. Resulting matrix values (mathematical representation)
3. Visual effects on 3D objects (rendered output)

## Overall Layout

```
┌─────────────────────────────────────────────────────────────┐
│                    Header (60px fixed)                       │
│  [☰ Matrix Demo]    [WebGL|Three.js]    [Theme] [← Home]   │
├────────────────┬────────────────────────────────────────────┤
│                │                                             │
│   Sidebar      │           Canvas Area                       │
│   200-900px    │           (auto-fill)                       │
│   Resizable    │                                             │
│                │                                             │
│  [Shapes]      │                                             │
│  [Controls]    │                                             │
│  [Panels...]   │                                             │
│                │                                             │
├────────────────┴────────────────────────────────────────────┤
│         Bottom Matrix Panel (200-600px resizable)           │
│   [Model] → [View] → [Projection] → [MVP] → [Screen]       │
└─────────────────────────────────────────────────────────────┘
```

## Component Specifications

### 1. Header (Topbar)
- **Height**: 60px (increased from 44px for breathing room)
- **Layout**: Flex row, space-between
- **Left**: Mobile menu button + Title "🔢 Matrix Demo"
- **Center**: Tab switcher (WebGL / Three.js)
- **Right**: Theme dropdown + "← 返回主页" link
- **Background**: `--bg-secondary` with bottom border
- **Sticky**: Fixed to top, z-index 50

#### Theme Switcher
- **Type**: Dropdown menu button
- **Position**: Header right side, before "Home" link
- **Style**: Same as other header buttons
- **Content**: 8 theme options (see Theme System section)
- **State**: Persisted to localStorage as `theme-preference`

#### Tab Switcher
- **Padding**: 8px 18px (per button)
- **Gap**: 2px between buttons
- **Border-radius**: 6px container, 4px buttons
- **Active state**:
  - Background: `--bg-tertiary`
  - Color: `--accent-primary`
  - Glow effect: `box-shadow: 0 0 12px var(--accent-primary-alpha)`
- **Mobile**: Hide text labels below 640px breakpoint, show icons only

### 2. Sidebar (Control Panel)
- **Width**: 200px - 900px (default 280px)
- **Background**: `--bg-secondary`
- **Border**: 1px solid `--border-color` on right
- **Overflow**: Vertical scroll, horizontal hidden
- **Padding**: 16px
- **State**: Width persisted to localStorage

#### Resize Handle
- **Width**: 4px
- **Position**: Absolute, at right edge of sidebar
- **Cursor**: col-resize
- **Default state**: Transparent background (subtle)
- **Hover/Dragging state**: `--accent-primary` background (prominent)
- **Transition**: background 0.15s ease

#### Shape Selector
- **Layout**: Horizontal button group
- **Button size**: Responsive based on sidebar width
  - Min: 40px × 40px
  - Max: 60px × 60px
  - Gap: 8px
- **Content**: SVG icon + text label below
- **Active state**:
  - Background: `--accent-primary` with 10% opacity
  - Border: 2px solid `--accent-primary`
  - Transform: translateY(-2px) (lift effect)
  - Box-shadow: 0 4px 12px `--accent-primary-alpha` (glow)
- **Shapes**: Cube, Sphere, Cone, Torus, Cylinder

#### Collapsible Panels
- **Spacing**: 16px between panels (increased from 12px)
- **Border**: 1px solid `--border-subtle`
- **Border-radius**: 6px
- **Background**: rgba(255, 255, 255, 0.01)

**Visual Hierarchy Modes** (user-toggleable):
1. **Importance mode**:
   - Primary panels (Model Transform): Thicker border, brighter background
   - Secondary panels: Standard styling
2. **Color-coded mode**:
   - Transform panels: Blue accent tint
   - Camera panels: Green accent tint
   - Rendering panels: Purple accent tint

**Summary (collapsed state)**:
- **Padding**: 8px 12px
- **Font**: 0.7rem, 600 weight, uppercase, 0.08em letter-spacing
- **Color**: `--text-secondary` (hover: `--text-primary`)
- **Chevron**: 12px SVG, rotates 180deg when open

**Content (expanded state)**:
- **Padding**: 10px 12px 12px
- **Animation**: slideDown 0.2s ease

**Default state**:
- Most important 1-2 panels open (e.g., "Model Transform")
- Others collapsed
- State persisted via localStorage

#### Slider Controls
- **Layout**:
  ```
  [Label left-aligned]              [Value right-aligned]
  [─────────●─────────────────────────────────────────]
  ```
- **Label**: 0.75rem, `--text-secondary`
- **Value**: 0.75rem, `--text-primary`, monospace
- **Track**:
  - Height: 6px
  - Background: `--bg-quaternary`
  - Border-radius: 3px
  - Filled portion: Linear gradient using `--accent-primary`
- **Thumb**:
  - 16px × 16px circle
  - Background: `--accent-primary`
  - Border: 2px solid `--bg-secondary`
  - Box-shadow: 0 2px 6px rgba(0, 0, 0, 0.3)
- **Spacing**: 8px between label row and slider, 16px between slider groups

### 3. Canvas Area
- **Size**: Auto-fill remaining space
- **Background**: `--bg-primary`
- **Content**: WebGL or Three.js renderer
- **Transition**: Smooth object transforms when parameters change

### 4. Bottom Matrix Panel
- **Height**: 200px - 600px (default 280px)
- **Width**: Full width (100%)
- **Background**: `--bg-secondary`
- **Border-top**: 1px solid `--border-color`
- **Overflow**: Horizontal scroll
- **State**: Height persisted to localStorage
- **Grid row**: 4 (fixed to bottom of layout)

#### Resize Handle
- **Height**: 8px (increased for better touch targets)
- **Position**: Absolute, at top edge of panel
- **Cursor**: ns-resize
- **Handle indicator**: 48px × 4px pill shape, centered
- **Default state**: `--border-color` background (subtle)
- **Hover/Dragging state**: `--accent-primary` background, expands to 64px × 5px
- **Glow effect**: `box-shadow: var(--shadow-glow)` when active
- **Transition**: all 0.15s ease
- **Mobile**: Hidden on devices ≤768px (prevents accidental touches)

#### Pipeline Flow Layout
Horizontal flowchart-style layout with arrows connecting each stage:

```
[Model Matrix] → [View Matrix] → [Projection Matrix] → [MVP Matrix]
```

**Matrix Card** (Final implementation):
- **Size**: 280px - 340px (desktop), 150px - 170px (mobile)
- **Border-radius**: 12px (`--radius-lg`)
- **Background**: `--bg-tertiary`
- **Border**: 1px solid `--border-color`
- **Padding**: `--space-xl` (24px desktop, 12px mobile)
- **Gap**: `--space-md` between internal elements
- **Box-shadow**: 0 2px 8px rgba(0, 0, 0, 0.2)
- **Hover**: Border becomes `--text-muted`, shadow increases

**Card Structure**:
```
┌─────────────────────────┐
│ Local → × Model → World │  ← Header (text-base/lg, colored stages)
│   T · Ry · Rx · Rz · S   │  ← Formula (text-md, background chip)
│                         │
│ ┌                     ┐ │  ← Math brackets (6rem)
│ │ 1.00  0.00  0.00  1 │ │  ← 4×4 grid, monospace
│ │ 0.00  1.00  0.00  0 │ │     1.125rem (desktop), 0.75rem (mobile)
│ │ 0.00  0.00  1.00  0 │ │     2 decimal places, right-aligned
│ │ 0.00  0.00  0.00  1 │ │     Tabular nums for alignment
│ └                     ┘ │
└─────────────────────────┘
```

**Matrix Card Header**:
- **From stage**: `--accent-secondary` color, `text-lg` size, bold weight
- **Operation**: `--text-muted` color, `text-md` size, monospace font
- **To stage**: `--accent-primary` color, `text-lg` size, bold weight
- **Spacing**: `--space-xs` between elements

**Matrix Formula**:
- **Font**: Monospace, `text-md` size
- **Color**: `--text-secondary`
- **Background**: `--bg-quaternary`
- **Padding**: 6px `--space-md` (4px on mobile)
- **Border-radius**: `--radius-sm`
- **Text-align**: Center

**Matrix Data Display**:
- **Font**: Cascadia Code / Consolas, 1.125rem (desktop), 0.75rem (mobile)
- **Grid**: 4×4 CSS Grid
  - Columns: `repeat(4, minmax(0, 1fr))` for equal widths
  - Rows: `repeat(4, auto)`
  - Gap: 4px 8px (row column), 2px 4px (mobile)
- **Precision**: Fixed 2 decimal places
- **Alignment**: Right-aligned, tabular nums
- **Positive numbers**: Leading space for alignment with negatives
- **Brackets**: 6rem font-size (desktop), 2.5rem (mobile)
- **Colors**:
  - Default: `--accent-primary`
  - Diagonal (highlighted): Same color, `font-weight: medium`
  - Background: `--bg-primary`

**Highlighting States**:
- **Default**: Standard border and background
- **Highlight active**:
  - Border: `--accent-primary`
  - Background: `--accent-primary-alpha` (10% opacity)
  - Box-shadow: `var(--shadow-glow-strong)`
  - Animation: `matrixPulse` (scale 1 → 1.03 → 1 over 0.6s)
- **Cell highlight**:
  - Animation: `cellHighlight` (scale 1 → 1.1 → 1 over 0.4s)
  - Used for diagonal elements and parameter correlations

**Arrow Connectors**:
- **Size**: 40px × 24px (desktop), 20px × 14px (mobile)
- **Color**: `--text-muted` (default)
- **Active state**: `--accent-primary`
- **Animation**: `arrowFlow` (1s ease-in-out infinite)
  - Opacity: 0.4 → 1 → 0.4
  - Transform: translateX(0) → translateX(4px) → translateX(0)

### 5. Pipeline Bar (Transform Visualization)
- **Height**: Auto (content-based), typically ~48px
- **Background**: #13191f (slightly darker than panel)
- **Border**: 1px solid `--border-color` (top and bottom)
- **Padding**: `--space-md` `--space-lg`
- **Overflow-x**: Auto
- **Layout**: Flex row, center-aligned
- **Gap**: `--space-md` between elements
- **Grid row**: 3

**Stage Labels (`.ps`)**:
- **Font-size**: 1rem (desktop), 0.85rem (mobile)
- **Font-weight**: 600
- **Padding**: 6px 16px (desktop), 4px 10px (mobile)
- **Background**: `--bg-secondary`
- **Border**: 1px solid `--border-color`
- **Border-radius**: `--radius-md`
- **Color**: `--text-secondary`
- **Small text**: 0.85rem (0.7rem mobile), `--accent-blue`

**Operation Labels (`.pa`)**:
- **Font-size**: 1rem (desktop), 0.85rem (mobile)
- **Padding**: 0 `--space-sm` (0 4px mobile)
- **Color**: `--text-secondary`
- **Strong**: 1.1rem (0.9rem mobile), `--accent-green`, 700 weight
- **Arrows**: 1.2rem (1rem mobile), `--text-muted`
- **Active state** (when parameter changes):
  - Particles flow from source matrix to destination
  - Color: `--accent-primary` with gradient fade
  - Duration: 0.8s
  - 8-12 particles per flow
- **Path**: Bezier curve or straight line depending on spacing

#### Highlight Behavior (Key Feature)

When user adjusts a parameter (e.g., "Rotation X" slider):

1. **Matrix Card Border**:
   - Affected matrix card (e.g., Model Matrix) gets glowing border
   - Color: `--accent-primary`
   - Animation: fade in 0.2s, pulse subtly, fade out 0.3s after adjustment stops

2. **Precise Value Highlighting**:
   - Specific matrix elements that changed get highlighted
   - Color: `--accent-primary`
   - Animation: Scale 1.0 → 1.2 → 1.0 over 0.4s (pulse effect)
   - Only affected values animate (e.g., rotation affects rotation matrix elements)

3. **Particle Flow**:
   - Particles flow from parameter control → affected matrix → canvas
   - Shows causality chain visually
   - Triggered on parameter change

4. **Canvas Transition**:
   - 3D object transforms smoothly (0.2s ease transition)

### 5. Buttons
**Primary Button** (e.g., "View Source"):
- **Padding**: 8px 16px
- **Background**: `--accent-primary` with 15% opacity
- **Border**: 1px solid `--accent-primary`
- **Color**: `--accent-primary`
- **Border-radius**: 6px
- **Hover**: Background opacity 25%, slight lift

**Secondary Button** (e.g., "Reset"):
- **Padding**: 8px 16px
- **Background**: `--bg-tertiary`
- **Border**: 1px solid `--border-color`
- **Color**: `--text-secondary`
- **Border-radius**: 6px
- **Hover**: Background `--bg-quaternary`, color `--text-primary`

## Interaction Patterns

### Parameter → Matrix → Visual Chain
This is the core educational feature. The design makes this relationship obvious through smart highlighting:

1. **User drags slider** (e.g., Rotation X: 0° → 45°)
2. **Immediate feedback**:
   - Slider value updates in real-time
   - **Matrix cards highlight** based on parameter type:
     - Model params (tx, ty, tz, rx, ry, rz, sx, sy, sz): All 4 matrices glow
     - Camera params (ex, ey, ez, ltx, lty, ltz): Last 3 matrices glow
     - Projection params (fov, near, far, ortho-size): Last 2 matrices glow
   - **Arrow connectors animate** showing data flow between matrices
   - **Matrix cards pulse** with scale animation (1 → 1.03 → 1 over 0.6s)
   - 3D object transforms smoothly
3. **Feedback timing**:
   - Matrix highlight: instant
   - Arrow animation: 1s ease-in-out infinite loop
   - Card pulse: 0.6s ease-out
   - Auto-clear: 2 seconds after last input
4. **Educational benefit**:
   - Shows causality: "changing this parameter affects these matrices"
   - Visualizes pipeline: data flows from Model → View → Projection → MVP
   - Helps students understand which transformations happen at each stage

### Collapsible Panel Behavior
- Click summary to toggle
- Chevron rotates 180° when opening
- Content slides down with opacity fade-in (slideDown animation, 0.2s)
- State saved to localStorage with unique key per panel
- Most important panels open by default (Model Transform, Camera)

### Resize Interactions
- **Sidebar resize**: Drag right edge, constrain to 200-900px
  - Modifies CSS grid column width
  - Visual feedback: Blue accent line during drag
  - Cursor: col-resize
- **Bottom panel resize**: Drag top edge, constrain to 200-600px
  - Modifies CSS grid row height
  - Handle expands and glows when hovering/dragging
  - Cursor: ns-resize
- Both save dimensions to localStorage
- Resize handles have dual-state visibility (subtle by default, prominent on hover)

### Theme Switching
- Click theme button in header right
- Dropdown with 8 theme options:
  - Color swatch preview
  - Theme name
  - "默认" badge on default theme
- Click theme to apply instantly
- Updates `data-theme` attribute on `<html>` root
- Preference saved to `localStorage: theme-preference`
- All CSS variables cascade automatically

## Responsive Breakpoints

### Desktop (> 1024px)
- Full three-panel layout (sidebar + canvas + matrix panel)
- Sidebar: 200-900px resizable
- Matrix cards: 280-340px wide
- All text labels visible
- Both resize handles active

### Tablet (769px - 1024px)
- Sidebar: 200-240px constrained
- Matrix cards: 200-240px wide
- Font sizes: ~90% of desktop
- All resize handles active
- Optimized spacing for medium screens

### Mobile (≤768px)
- Single column layout
- Sidebar becomes fixed overlay (280px wide, slides from left)
- Canvas area: minmax(300px, 1fr)
- Pipeline bar: auto height
- Matrix panel: 200px fixed (not resizable on mobile)
- Matrix cards: 150-170px wide
- Font sizes: ~65% of desktop (0.75rem for matrix values)
- Bottom panel resize handle hidden (prevents accidental touches)
- Touch-friendly targets: min 44px

### Small Mobile (≤480px)
- Matrix pipeline: wraps to center, allows scrolling
- Matrix cards: 130-150px wide
- Font sizes: Further reduced (0.7rem for matrix values)
- Brackets: 2rem (smaller for compact display)
- Optimized for very narrow screens

### Touch Optimizations
- Larger touch targets (min 44px × 44px)
- No hover-dependent interactions
- Simplified animations (reduced motion support)
- Scrollable areas with momentum scrolling

## Accessibility
- All interactive elements have ARIA labels
- Keyboard navigation supported:
  - Tab: Navigate between controls
  - Enter/Space: Activate buttons and toggles
  - Escape: Close dropdowns and overlays
- Focus indicators: 2px solid `--accent-primary`, 1px offset
- Color contrast: WCAG AA compliant (4.5:1 for normal text)
- Resize handles: `role="separator"`, `aria-orientation`
- Theme picker: `role="menu"`, `aria-haspopup`, `aria-expanded`
- Dropdown menus: Keyboard accessible (arrow keys to navigate)
- Reduced motion support: Respects `prefers-reduced-motion`

## Performance Considerations
- CSS transitions: Hardware-accelerated (transform, opacity)
- Animations: GPU-accelerated properties only
- Matrix calculations: Throttled to 60fps max
- LocalStorage writes: Immediate for critical settings
- Highlight timeouts: Debounced (2s auto-clear)
- Responsive images: Future consideration for matrix visualizations
- Code splitting: Astro handles route-based splitting
- Canvas rendering optimized (only redraw on change)
