import urllib.request
import os
import math
import numpy as np

# -------------------------
# Global definitions
# -------------------------
aromatic_aa = {"HIS", "TYR", "TRP", "PHE"}
AROMATIC_RINGS = {
    "HIS": [{"CG", "ND1", "CD2", "CE1", "NE2"}],       # imidazole ring
    "PHE": [{"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}], # phenyl ring
    "TYR": [{"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}], # phenol ring
    "TRP": [{"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"}] # benzene ring of TRP
}


# -------------------------
# Utility Functions
# -------------------------
def download_pdb_urllib(pdb_id, save_dir="."):
    """
    Download a PDB file using urllib.
    Returns the local file path.
    """
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb"
    file_path = os.path.join(save_dir, f"{pdb_id.lower()}.pdb")

    try:
        urllib.request.urlretrieve(url, file_path)
        print(f"Successfully downloaded: {file_path}")
        return file_path
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        return None


def extract_atom_positions(pdb_id_file):
    """
    Extract atom coordinates from a PDB file.
    Returns a list of tuples:
    (residue_name, residue_seq, atom_name, x, y, z)
    Only includes atoms in chain A.
    """
    positions = []

    with open(pdb_id_file, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    residue_seq = int(line[22:26].strip())
                    chain_id = line[21:22].strip()
                    if chain_id == "A":
                        positions.append(
                            (residue_name, residue_seq, atom_name, x, y, z)
                        )
                except ValueError:
                    continue
    return positions


def plane_from_points(points):
    """
    Given a list of (x, y, z) coordinates, compute the plane coefficients
    a, b, c, d such that ax + by + cz + d = 0.
    Returns (a, b, c, d, centroid)
    """
    pts = np.array(points)

    # Compute centroid
    centroid = pts.mean(axis=0)

    # Shift points relative to centroid
    shifted = pts - centroid

    # Compute normal vector using SVD
    _, _, vh = np.linalg.svd(shifted)
    normal = vh[-1]

    a, b, c = normal
    d = -np.dot(normal, centroid)

    return a, b, c, d, centroid


def compute_aromatic_coordinates(aromatic_atoms):
    """
    Compute the centroid and plane of aromatic rings.
    Returns a list of tuples:
    (resname, resnum, cx, cy, cz, a, b, c, d)
    """
    from collections import defaultdict

    residues = defaultdict(list)

    # Group atoms by residue
    for resname, resnum, atomname, x, y, z in aromatic_atoms:
        if resname in AROMATIC_RINGS:
            residues[(resname, resnum)].append((atomname, x, y, z))

    aromatic_list = []

    for (resname, resnum), atoms in residues.items():
        atom_dict = {a: (x, y, z) for a, x, y, z in atoms}
        ring_atoms = AROMATIC_RINGS[resname][0]

        # Ensure complete set of atoms exists
        if not all(a in atom_dict for a in ring_atoms):
            continue

        coords = [atom_dict[a] for a in ring_atoms]

        # Compute plane and centroid
        a, b, c, d, centroid = plane_from_points(coords)
        cx, cy, cz = centroid

        aromatic_list.append(
            (resname, resnum, cx, cy, cz, a, b, c, d)
        )

    return aromatic_list


def angle_between_line_and_plane(cys_coordinates, aromatic_coordinates,
                                cutoff_outer=10.5, cutoff_inner=3.0):
    """
    Compute distances and angles between cysteine sulfur atoms and
    aromatic ring centers. Returns a list with detailed info including
    plane coefficients and line vectors.
    """
    results = []

    # Header describing each field
    results.append((
        "cys_resname",
        "cys_resnum",
        "aro_resname",
        "aro_resnum",
        "dist from S to aro center",
        "angle of line and plane",
        "Lx", "Ly", "Lz",
        "sx", "sy", "sz",
        "plane_a", "plane_b", "plane_c", "plane_d"
    ))

    cutoff_outer_sq = cutoff_outer ** 2
    cutoff_inner_sq = cutoff_inner ** 2

    for cys in cys_coordinates:
        cys_resname = cys[0]
        cys_resnum = cys[1]
        sx, sy, sz = cys[3], cys[4], cys[5]

        S = np.array([sx, sy, sz])

        for aro in aromatic_coordinates:
            aro_resname, aro_resnum = aro[0], aro[1]
            cx, cy, cz = aro[2], aro[3], aro[4]
            a, b, c, d = aro[5], aro[6], aro[7], aro[8]

            C = np.array([cx, cy, cz])
            L = C - S

            # Distance filter
            dist_sq = np.sum(L**2)
            if dist_sq > cutoff_outer_sq or dist_sq < cutoff_inner_sq:
                continue

            # Angle calculation
            N = np.array([a, b, c])
            L_norm = np.linalg.norm(L)
            N_norm = np.linalg.norm(N)
            if L_norm == 0 or N_norm == 0:
                continue

            dot_val = abs(np.dot(L, N))
            phi = math.degrees(math.acos(dot_val / (L_norm * N_norm)))
            theta = 90.0 - phi

            # Apply rounding
            distance = round(math.sqrt(dist_sq), 1)
            theta = round(theta)

            # Line vector
            Lx, Ly, Lz = L

            # Append full info
            results.append((
                cys_resname, cys_resnum,
                aro_resname, aro_resnum,
                distance, theta,
                Lx, Ly, Lz,
                sx, sy, sz,
                a, b, c, d
            ))

    return results


def point_to_trimmed_segment_distance(P, A, B, trim=0.6):
    """
    Compute distance from point P to a line segment AB, trimmed
    0.6 Å from each end to avoid endpoints.
    """
    Px, Py, Pz = P
    Ax, Ay, Az = A
    Bx, By, Bz = B

    AB = (Bx - Ax, By - Ay, Bz - Az)
    AB_len = math.sqrt(sum([i**2 for i in AB]))
    if AB_len <= trim / 2:
        return float("inf")

    ux, uy, uz = (AB[0] / AB_len, AB[1] / AB_len, AB[2] / AB_len)
    A2 = (Ax + trim * ux, Ay + trim * uy, Az + trim * uz)
    B2 = (Bx - trim * ux, By - trim * uy, Bz - trim * uz)

    AB2 = (B2[0] - A2[0], B2[1] - A2[1], B2[2] - A2[2])
    AB2_len2 = sum([i**2 for i in AB2])
    AP = (Px - A2[0], Py - A2[1], Pz - A2[2])

    t = sum(AP[i] * AB2[i] for i in range(3)) / AB2_len2
    t_clamped = max(0, min(1, t))

    C = tuple(A2[i] + t_clamped * AB2[i] for i in range(3))
    return math.sqrt(sum((P[i] - C[i])**2 for i in range(3)))


def determine_obstruction(results, all_atoms, threshold=0.5):
    """
    Determine if any atom lies within threshold of the S→aromatic center line segment.
    Returns updated results with 'obstructed' or 'un-obstructed' appended.
    """
    updated_results = []
    header = list(results[0]) + ["obstruction"]
    updated_results.append(tuple(header))

    for entry in results[1:]:
        (cys_resname, cys_resnum,
         aro_resname, aro_resnum,
         distance, angle,
         Lx, Ly, Lz,
         sx, sy, sz,
         a, b, c, d) = entry

        A = (sx, sy, sz)
        B = (Lx, Ly, Lz)
        obstructed = False

        for atom in all_atoms:
            _, _, _, ax, ay, az = atom
            P = (ax, ay, az)
            d_seg = point_to_trimmed_segment_distance(P, A, B)
            if d_seg <= threshold:
                obstructed = True
                break

        updated_results.append(entry + ("obstructed" if obstructed else "un-obstructed",))

    return updated_results


def categorize_angle(angle):
    """Convert numeric angle into sigma/mid/pi string."""
    if 0 <= angle < 20:
        return "sigma"
    elif 20 <= angle < 70:
        return "mid"
    elif 70 <= angle <= 90:
        return "pi"
    return "out-of-range"


def refine_final_list(results):
    """
    Refine the results:
    - remove coordinate and plane data
    - convert numeric angles to sigma/mid/pi
    Returns final clean list.
    """
    refined = []
    refined.append((
        "cys_resname",
        "cys_resnum",
        "aro_resname",
        "aro_resnum",
        "distance",
        "angle_type",
        "obstruction"
    ))

    for entry in results[1:]:
        (cys_resname, cys_resnum,
         aro_resname, aro_resnum,
         distance, angle,
         Lx, Ly, Lz,
         sx, sy, sz,
         a, b, c, d,
         obstruction) = entry

        angle_type = categorize_angle(angle)

        refined.append((
            cys_resname,
            cys_resnum,
            aro_resname,
            aro_resnum,
            distance,
            angle_type,
            obstruction
        ))

    return refined


# -------------------------
# Main Pipeline Function
# -------------------------
def run_protein_structure_solver(pdb_id):
    """
    Run full analysis pipeline:
    1. Download PDB
    2. Extract atoms
    3. Compute aromatic centers
    4. Compute distances and angles
    5. Check linear obstruction
    6. Refine final output
    Returns final clean list.
    """
    pdb_id_file = download_pdb_urllib(pdb_id)
    all_atoms = extract_atom_positions(pdb_id_file)

    aromatic_atoms = [atom for atom in all_atoms if atom[0] in aromatic_aa]
    cys_atoms = [atom for atom in all_atoms if atom[0] == "CYS"]
    cys_coordinates = [atom for atom in cys_atoms if atom[2] == "SG"]

    aromatic_coordinates = compute_aromatic_coordinates(aromatic_atoms)

    full_distance_and_angle = angle_between_line_and_plane(
        cys_coordinates, aromatic_coordinates
    )

    aromatic_data_with_obstruction = determine_obstruction(
        full_distance_and_angle, all_atoms
    )

    final_clean_list = refine_final_list(aromatic_data_with_obstruction)
    return final_clean_list


# -------------------------
# Execute analysis
# -------------------------
final_output = run_protein_structure_solver("4F5S")
