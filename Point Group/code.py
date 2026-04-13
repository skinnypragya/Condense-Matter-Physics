import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle

# --- Style matching the LaTeX document ---
SHADE_BLUE = '#5340B7'
MIRROR_ORANGE = '#EF9F27'
ROT_BLUE = '#3B8BD4'
SHAPE_FILL = '#F0F0F0'
LINE_COLOR = '#999999'

# --- Geometry Generators ---
def get_square_slices():
    """Generates the 8 triangular slices of the base square, clockwise from top."""
    slices = []
    # Vertices clockwise starting from top-center
    pts = [(0,1), (1,1), (1,0), (1,-1), (0,-1), (-1,-1), (-1,0), (-1,1)]
    for i in range(8):
        slices.append(np.array([[0, 0], pts[i], pts[(i+1)%8]]))
    return slices

def get_hex_slices():
    """Generates the 12 triangular slices of the base hexagon, clockwise from top."""
    slices = []
    # Pointy-top hexagon vertices
    v = [(np.cos(a), np.sin(a)) for a in np.radians([90, 30, -30, -90, -150, 150])]
    pts = []
    for i in range(6):
        pts.append(v[i]) # Vertex
        # Edge Midpoint
        m = ((v[i][0] + v[(i+1)%6][0])/2, (v[i][1] + v[(i+1)%6][1])/2)
        pts.append(m)
    for i in range(12):
        slices.append(np.array([[0, 0], pts[i], pts[(i+1)%12]]))
    return slices

# --- Point Group Definitions ---
# "shaded": indices of slices to shade (0-indexed). 
# "ops": List of symmetry operations (type, angle_in_degrees)
groups = {
    "C1":  {"base": "square", "shaded": [0], "marker": None, "mirrors": [],
            "ops": [("rot", 360)]},
    "C2":  {"base": "square", "shaded": [0, 4], "marker": "circle", "mirrors": [],
            "ops": [("rot", 180), ("rot", 360)]},
    "Cs":  {"base": "square", "shaded": [0, 7], "marker": None, "mirrors": [90],
            "ops": [("rot", 360), ("ref", 90)]},
    "C2v": {"base": "square", "shaded": [0, 3, 4, 7], "marker": "circle", "mirrors": [0, 90],
            "ops": [("rot", 180), ("rot", 360), ("ref", 0), ("ref", 90)]},
    "C4":  {"base": "square", "shaded": [0, 2, 4, 6], "marker": "square", "mirrors": [],
            "ops": [("rot", 90), ("rot", 180), ("rot", 270), ("rot", 360)]},
    "C4v": {"base": "square", "shaded": [], "marker": "square", "mirrors": [0, 45, 90, 135],
            "ops": [("rot", 90), ("rot", 180), ("rot", 270), ("rot", 360),
                    ("ref", 0), ("ref", 45), ("ref", 90), ("ref", 135)]},
    
    "C3":  {"base": "hex", "shaded": [0, 4, 8], "marker": "triangle", "mirrors": [],
            "ops": [("rot", 120), ("rot", 240), ("rot", 360)]},
    "C3v": {"base": "hex", "shaded": [0, 11, 3, 4, 7, 8], "marker": "triangle", "mirrors": [30, 90, 150],
            "ops": [("rot", 120), ("rot", 240), ("rot", 360), ("ref", 30), ("ref", 90), ("ref", 150)]},
    "C6":  {"base": "hex", "shaded": [0, 2, 4, 6, 8, 10], "marker": "hexagon", "mirrors": [],
            "ops": [("rot", 60), ("rot", 120), ("rot", 180), ("rot", 240), ("rot", 300), ("rot", 360)]},
    "C6v": {"base": "hex", "shaded": [], "marker": "hexagon", "mirrors": [0, 30, 60, 90, 120, 150],
            "ops": [("rot", 60*i) for i in range(1, 7)] + [("ref", 30*i) for i in range(6)]}
}

# --- Mathematical Symmetry Check Engine ---
def rotate_point(p, angle_deg):
    th = np.radians(angle_deg)
    c, s = np.cos(th), np.sin(th)
    return np.array([p[0]*c - p[1]*s, p[0]*s + p[1]*c])

def reflect_point(p, angle_deg):
    th = np.radians(angle_deg)
    c, s = np.cos(th), np.sin(th)
    u = np.array([c, s])
    proj = np.dot(p, u)
    return 2 * proj * u - p

def perform_symmetry_operations():
    """Mathematically proves the symmetry by applying ops to the centroids of shaded regions."""
    sq_slices = get_square_slices()
    hx_slices = get_hex_slices()
    
    print("=" * 60)
    print("PERFORMING AND VERIFYING SYMMETRY OPERATIONS")
    print("=" * 60)
    
    for name, config in groups.items():
        slices = sq_slices if config["base"] == "square" else hx_slices
        shaded_idx = config["shaded"] if config["shaded"] else list(range(len(slices)))
        
        # Calculate centroids of all shaded sectors
        centroids = np.array([np.mean(slices[i], axis=0) for i in shaded_idx])
        
        print(f"\n[{name}] Group -> Testing {len(config['ops'])} Operations:")
        for op_type, val in config["ops"]:
            valid = True
            for c in centroids:
                # Apply transformation matrix
                c_trans = rotate_point(c, val) if op_type == "rot" else reflect_point(c, val)
                
                # Verify the transformed point perfectly overlaps an existing shaded centroid
                distances = np.linalg.norm(centroids - c_trans, axis=1)
                if np.min(distances) > 1e-4:
                    valid = False
                    break
                    
            op_str = f"{'Rotation' if op_type=='rot' else 'Reflection'} @ {val}°"
            # print(f"   |-- {op_str:<18} : {'[✓] INVARIANT' if valid else '[✗] BROKEN'}")

# --- Drawing Engine ---
def draw_shape(ax, slices, shaded_idx, title, marker, mirrors):
    # 1. Draw geometric slices
    for i, poly in enumerate(slices):
        color = SHADE_BLUE if i in shaded_idx else SHAPE_FILL
        ax.add_patch(Polygon(poly, facecolor=color, edgecolor=LINE_COLOR, linewidth=0.6, zorder=1))
    
    # 2. Draw mirror planes
    for ang in mirrors:
        rad = np.radians(ang)
        dx, dy = np.cos(rad), np.sin(rad)
        ax.plot([-1.1*dx, 1.1*dx], [-1.1*dy, 1.1*dy], color=MIRROR_ORANGE, linestyle='--', lw=1.5, zorder=2)
        
    # 3. Draw rotational markers
    if marker == 'circle':
        ax.add_patch(Circle((0, 0), 0.08, facecolor=ROT_BLUE, edgecolor='white', lw=1, zorder=3))
    elif marker == 'square':
        ax.add_patch(Polygon([[-0.08, -0.08], [0.08, -0.08], [0.08, 0.08], [-0.08, 0.08]], 
                             facecolor=ROT_BLUE, edgecolor='white', lw=1, zorder=3))
    elif marker == 'triangle':
        ax.add_patch(Polygon([[0, 0.1], [-0.09, -0.05], [0.09, -0.05]], 
                             facecolor=ROT_BLUE, edgecolor='white', lw=1, zorder=3))
    elif marker == 'hexagon':
        pts = [[0.1*np.cos(a), 0.1*np.sin(a)] for a in np.radians([90, 30, -30, -90, -150, 150])]
        ax.add_patch(Polygon(pts, facecolor=ROT_BLUE, edgecolor='white', lw=1, zorder=3))

    # Format canvas
    ax.set_aspect('equal')
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.axis('off')
    ax.set_title(title, y=-0.1, fontsize=12)

def plot_all_groups():
    sq_slices = get_square_slices()
    hx_slices = get_hex_slices()
    
    fig, axes = plt.subplots(2, 5, figsize=(16, 7))
    axes = axes.flatten()
    
    # Titles mapping for presentation
    titles = {"C1":"$C_1$ (1)", "C2":"$C_2$ (2)", "Cs":"$C_s$ (m)", "C2v":"$C_{2v}$ (2mm)", 
              "C4":"$C_4$ (4)", "C4v":"$C_{4v}$ (4mm)", "C3":"$C_3$ (3)", "C3v":"$C_{3v}$ (3m)", 
              "C6":"$C_6$ (6)", "C6v":"$C_{6v}$ (6mm)"}
    
    for ax, (name, config) in zip(axes, groups.items()):
        slices = sq_slices if config["base"] == "square" else hx_slices
        draw_shape(ax, slices, config["shaded"], titles[name], config["marker"], config["mirrors"])

    plt.suptitle("2D Crystallographic Point Groups (Schoenflies & Hermann-Mauguin)", fontsize=16, y=0.98)
    plt.tight_layout()
    plt.show()

# --- Execution ---
if __name__ == "__main__":
    # 1. Output the mathematical execution of the operators
    perform_symmetry_operations()
    
    # 2. Render the graphics 
    plot_all_groups()