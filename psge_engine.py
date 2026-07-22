"""
psge_engine.py — PSGE II v2.6.1 : Universe UX Stabilization
============================================================
Discrete Simplicial Geometry Engine for Regge Calculus.

Public API (stable since v2.6.1):
  - MetricTensor      : generalised metric tensor (Euclidean / Lorentzian)
  - MetricType        : EUCLIDEAN | LORENTZIAN
  - GeometryEngine    : low-level orchestrator (certificate + coordinates → report)
  - Universe          : high-level UX wrapper for a complete discrete universe
  - ReggeGeometryReport : immutable result with O(1) edge lookup
  - HingeTopology     : INTERIOR | BOUNDARY | NON_MANIFOLD | DEGENERATE
"""

__version__ = "2.6.1"

from dataclasses import dataclass, field
from typing import Optional, Dict, Tuple, List, Mapping
from types import MappingProxyType
from enum import Enum
from math import factorial
import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# 1. Exception Hierarchy
# ─────────────────────────────────────────────────────────────────────────────

class PSGEError(Exception):
    """Base class for all PSGE II specific errors."""


class TopologyError(PSGEError):
    """Raised when the topological certificate is invalid or inconsistent."""


class MetricError(PSGEError):
    """Raised when the metric tensor is degenerate or has an invalid signature."""


class GeometryError(PSGEError):
    """Raised when a geometric computation fails (degeneracy, insufficient rank)."""


# ─────────────────────────────────────────────────────────────────────────────
# 2. Enumerations
# ─────────────────────────────────────────────────────────────────────────────

class GeometryStatus(Enum):
    INVALID = 0
    COMPUTED = 1
    FAILED = 2


class MetricType(Enum):
    EUCLIDEAN = "Euclidean"
    LORENTZIAN = "Lorentzian"


class OrientationSign(Enum):
    POSITIVE = 1
    NEGATIVE = -1

    @classmethod
    def from_value(cls, val: float) -> "OrientationSign":
        return cls.POSITIVE if val >= 0 else cls.NEGATIVE


class GramSignature(Enum):
    EUCLIDEAN = "Euclidean"
    SPACELIKE = "Spacelike"
    TIMELIKE = "Timelike"
    DEGENERATE = "Degenerate"


class HingeTopology(Enum):
    """Topological classification of a hinge (edge) in the triangulation.

    INTERIOR     — the link of the edge is a closed 1-cycle (S¹); the edge is
                   surrounded by a complete ring of tetrahedra (bulk geometry).
    BOUNDARY     — the link is an open path; the edge sits on the boundary of
                   the manifold.
    NON_MANIFOLD — the link is neither a cycle nor a simple arc (degenerate
                   topology; multiple components, figure-eight, etc.).
    DEGENERATE   — no incident cells, or the link graph is empty.
    """
    INTERIOR = 1
    BOUNDARY = 2
    NON_MANIFOLD = 3
    DEGENERATE = 4


# ─────────────────────────────────────────────────────────────────────────────
# 3. Metric Tensor (generalised, immutable)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class MetricTensor:
    """Immutable, copy-defensive generalised metric tensor g_{ij}.

    All NumPy arrays stored inside are made read-only so that the frozen
    dataclass guarantee is honoured end-to-end.
    """
    _matrix: np.ndarray

    def __post_init__(self) -> None:
        matrix = np.array(object.__getattribute__(self, "_matrix"), dtype=float, copy=True)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise MetricError("The metric matrix must be square.")
        if not np.allclose(matrix, matrix.T):
            raise MetricError("The metric tensor must be symmetric.")
        matrix.setflags(write=False)
        object.__setattr__(self, "_matrix", matrix)

    # ── scalar products ───────────────────────────────────────────────────────

    def inner(self, u: np.ndarray, v: np.ndarray) -> float:
        """Generalised inner product g(u, v) = uᵀ G v."""
        return float(u.T @ self._matrix @ v)

    def norm2(self, v: np.ndarray) -> float:
        """Squared norm g(v, v)."""
        return self.inner(v, v)

    def norm(self, v: np.ndarray, eps: float = 1e-12) -> float:
        """Positive norm √max(0, g(v,v)), raises GeometryError if negative."""
        val = self.norm2(v)
        if val < -eps:
            raise GeometryError(
                f"Negative squared norm ({val:.6e}) incompatible with a Euclidean metric."
            )
        return float(np.sqrt(max(0.0, val)))

    def interval_sq(self, u: np.ndarray) -> float:
        """Spacetime interval squared g(u, u) (alias for norm2, Lorentzian use)."""
        return self.norm2(u)

    # ── derived operations ────────────────────────────────────────────────────

    def project_orthogonal(
        self, v: np.ndarray, direction: np.ndarray, eps: float = 1e-12
    ) -> np.ndarray:
        """Project *v* onto the hyperplane orthogonal to *direction*."""
        d2 = self.norm2(direction)
        if abs(d2) < eps:
            return np.zeros_like(v)
        return v - (self.inner(v, direction) / d2) * direction

    def gram_matrix(self, vectors: List[np.ndarray]) -> np.ndarray:
        """Gram matrix G_{ij} = g(vᵢ, vⱼ) for a list of vectors."""
        n = len(vectors)
        G = np.empty((n, n), dtype=float)
        for i in range(n):
            for j in range(n):
                G[i, j] = self.inner(vectors[i], vectors[j])
        return G

    def gram_area(self, u: np.ndarray, v: np.ndarray, eps: float = 1e-12) -> float:
        """Area of the parallelogram spanned by *u* and *v* via the metric."""
        g11 = self.inner(u, u)
        g22 = self.inner(v, v)
        g12 = self.inner(u, v)
        gram_det = g11 * g22 - g12 * g12
        if gram_det < -eps:
            raise GeometryError(f"Indefinite surface Gram determinant ({gram_det:.6e}).")
        return 0.5 * float(np.sqrt(max(0.0, gram_det)))

    # ── constructors ──────────────────────────────────────────────────────────

    @classmethod
    def euclidean(cls, dim: int) -> "MetricTensor":
        """Standard Euclidean metric δ_{ij} in *dim* dimensions."""
        return cls(np.eye(dim))

    @classmethod
    def minkowski(cls, dim: int = 4) -> "MetricTensor":
        """Minkowski metric diag(-1, 1, …, 1) in *dim* dimensions."""
        mat = np.eye(dim)
        mat[0, 0] = -1.0
        return cls(mat)


# ─────────────────────────────────────────────────────────────────────────────
# 4. Immutable Data Structures
# ─────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class IncidentCell:
    simplex: Tuple[int, ...]
    orientation: OrientationSign
    incident_faces: Tuple[Tuple[int, int, int], ...]
    opposite_vertices: Tuple[int, ...]


@dataclass(frozen=True)
class SimplexVolumeResult:
    volume: float
    gram_determinant: float
    signature: GramSignature
    is_valid: bool


@dataclass(frozen=True)
class HingeTraversalResult:
    hinge: Tuple[int, int]
    incident_cells: Tuple[IncidentCell, ...]
    dual_measure: float
    topology: HingeTopology


@dataclass(frozen=True)
class CurvatureReport:
    """Curvature data for a single hinge (edge) of the triangulation."""
    edge: Tuple[int, int]
    incident_cells: Tuple[IncidentCell, ...]
    deficit_angle: float
    hinge_measure: float
    dual_volume: float
    scalar_curvature: float
    topology: HingeTopology


@dataclass(frozen=True)
class GeometryStatistics:
    total_volume: float
    edge_lengths: Tuple[float, ...]
    gram_volumes: Tuple[float, ...]
    num_simplices: int
    num_degenerate_simplices_ignored: int


# ─────────────────────────────────────────────────────────────────────────────
# 5. Geometry Report (immutable, O(1) edge lookup)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class ReggeGeometryReport:
    """Complete Regge geometry report for a triangulated manifold.

    Attributes
    ----------
    bulk_regge_sum       : Regge action summed over interior hinges only.
    boundary_hinge_sum   : Regge action summed over boundary hinges only.
    total_geometric_sum  : Total Regge action (bulk + boundary).
    total_action         : Alias for *total_geometric_sum* (backward compatibility).
    """
    status: GeometryStatus
    total_volume: float
    bulk_regge_sum: float
    boundary_hinge_sum: float
    total_geometric_sum: float
    average_scalar_curvature: float
    maximum_deficit: float
    minimum_deficit: float
    total_deficit: float
    max_edge_length: float
    min_edge_length: float
    num_hinges_analyzed: int
    num_simplices: int
    num_degenerate_simplices_ignored: int
    min_gram_volume: float
    max_gram_volume: float
    curvature_reports: Tuple[CurvatureReport, ...]
    _edge_index_map: Mapping[Tuple[int, int], CurvatureReport] = field(
        default_factory=dict, init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        mapping = {report.edge: report for report in self.curvature_reports}
        object.__setattr__(self, "_edge_index_map", MappingProxyType(mapping))

    @property
    def total_action(self) -> float:
        """Backward-compatible alias for *total_geometric_sum*."""
        return self.total_geometric_sum

    def get_report(self, edge: Tuple[int, int]) -> Optional[CurvatureReport]:
        """O(1) lookup of a CurvatureReport by edge.  Edge vertices are sorted
        automatically so both ``(a, b)`` and ``(b, a)`` work."""
        return self._edge_index_map.get(tuple(sorted(edge)))

    def summary(self) -> str:
        """Human-readable summary string."""
        return (
            f"=== REGGE DISCRETE GEOMETRY REPORT (PSGE II v{__version__}) ===\n"
            f" Status                    : {self.status.name}\n"
            f" Total Volume              : {self.total_volume:.6f}\n"
            f" Regge Action (total)      : {self.total_geometric_sum:.6f}\n"
            f"   ↳ Bulk (interior hinges): {self.bulk_regge_sum:.6f}\n"
            f"   ↳ Boundary hinges       : {self.boundary_hinge_sum:.6f}\n"
            f" Average Scalar Curvature  : {self.average_scalar_curvature:.6f}\n"
            f" Deficit (Min/Max/Total)   : {self.minimum_deficit:.4f} / "
            f"{self.maximum_deficit:.4f} / {self.total_deficit:.4f}\n"
            f" Edge Length (Min/Max)     : {self.min_edge_length:.4f} / "
            f"{self.max_edge_length:.4f}\n"
            f" Hinges Analyzed           : {self.num_hinges_analyzed}\n"
            f" Simplices / Ignored       : {self.num_simplices} / "
            f"{self.num_degenerate_simplices_ignored}\n"
            f" Gram Volume (Min/Max)     : {self.min_gram_volume:.6f} / "
            f"{self.max_gram_volume:.6f}\n"
            f"{'=' * 52}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# 6. Topology Analysis
# ─────────────────────────────────────────────────────────────────────────────

class HingeTopologyAnalyzer:
    """Classify the topological type of a hinge from its incident cells.

    The *link* of a 1-simplex (edge) in a 3-manifold triangulation is a
    1-complex.  We distinguish:
      - S¹ (closed cycle, all vertices degree-2, Euler χ = 0) → INTERIOR
      - open arc (not all degree-2, but connected) → BOUNDARY
      - everything else → NON_MANIFOLD
    """

    @staticmethod
    def classify_link(incident_cells: Tuple[IncidentCell, ...]) -> HingeTopology:
        if not incident_cells:
            return HingeTopology.DEGENERATE

        link_edges: set = set()
        for cell in incident_cells:
            others = tuple(sorted(cell.opposite_vertices))
            if len(others) == 2:
                link_edges.add(others)
            elif len(others) > 2:
                for i in range(len(others)):
                    for j in range(i + 1, len(others)):
                        link_edges.add(tuple(sorted((others[i], others[j]))))

        degree: Dict[int, int] = {}
        for u, v in link_edges:
            degree[u] = degree.get(u, 0) + 1
            degree[v] = degree.get(v, 0) + 1

        if not degree:
            return HingeTopology.DEGENERATE

        # Build adjacency for connectivity check
        adjacency: Dict[int, List[int]] = {}
        for u, v in link_edges:
            adjacency.setdefault(u, []).append(v)
            adjacency.setdefault(v, []).append(u)

        # BFS connectivity
        start = next(iter(degree))
        visited: set = set()
        queue = [start]
        while queue:
            node = queue.pop(0)
            if node not in visited:
                visited.add(node)
                queue.extend(n for n in adjacency.get(node, []) if n not in visited)

        num_vertices = len(degree)
        num_edges = len(link_edges)
        euler_char = num_vertices - num_edges  # χ = 0 for S¹
        is_connected = len(visited) == num_vertices
        is_regular_deg2 = all(d == 2 for d in degree.values())
        is_cycle = num_edges == num_vertices

        if is_connected and is_regular_deg2 and is_cycle and euler_char == 0:
            return HingeTopology.INTERIOR
        elif is_connected and not is_regular_deg2:
            return HingeTopology.BOUNDARY
        else:
            return HingeTopology.NON_MANIFOLD


# ─────────────────────────────────────────────────────────────────────────────
# 7. Geometric Operators
# ─────────────────────────────────────────────────────────────────────────────

class SimplexVolume:
    """Compute the n-volume of a simplex via its Gram matrix."""

    @staticmethod
    def compute_volume(
        coords: Tuple[np.ndarray, ...],
        tensor: MetricTensor,
        metric_type: MetricType = MetricType.EUCLIDEAN,
        strict: bool = True,
    ) -> SimplexVolumeResult:
        v0 = coords[0]
        vectors = [v - v0 for v in coords[1:]]
        dim = len(vectors)
        gram = tensor.gram_matrix(vectors)

        rank = np.linalg.matrix_rank(gram)
        if rank < dim:
            if strict:
                raise GeometryError(
                    f"Degenerate simplex: Gram rank ({rank}) < dimension ({dim})."
                )
            return SimplexVolumeResult(0.0, 0.0, GramSignature.DEGENERATE, False)

        det = float(np.linalg.det(gram))

        if metric_type == MetricType.LORENTZIAN:
            sig = (
                GramSignature.SPACELIKE
                if det > 0
                else GramSignature.TIMELIKE
                if det < 0
                else GramSignature.DEGENERATE
            )
        else:
            if det <= 0:
                if strict:
                    raise GeometryError(
                        f"Non-positive Gram determinant ({det:.6e}) in Euclidean metric."
                    )
                return SimplexVolumeResult(0.0, det, GramSignature.DEGENERATE, False)
            sig = GramSignature.EUCLIDEAN

        volume = float(np.sqrt(abs(det)) / factorial(dim))
        return SimplexVolumeResult(volume, det, sig, True)


class BarycentricDualCell:
    """Barycentric dual measure for a hinge (edge).

    Sums the areas of the dual triangles formed at the edge midpoint,
    face barycenters and tetrahedral barycenters.
    """

    def compute_measure(
        self,
        edge: Tuple[int, int],
        incident_cells: Tuple[IncidentCell, ...],
        coords_provider,
        tensor: MetricTensor,
    ) -> float:
        a, b = edge
        pa = coords_provider.vertex(a)
        pb = coords_provider.vertex(b)

        if tensor.norm2(pb - pa) < 1e-15:
            return 0.0

        mid_ab = 0.5 * (pa + pb)
        dual_area_sum = 0.0
        for cell in incident_cells:
            tet_coords = [coords_provider.vertex(v) for v in cell.simplex]
            tet_bary = sum(tet_coords) / len(tet_coords)
            others = [v for v in cell.simplex if v != a and v != b]
            for c in others:
                pc = coords_provider.vertex(c)
                face_bary = (pa + pb + pc) / 3.0
                v1 = face_bary - mid_ab
                v2 = tet_bary - mid_ab
                dual_area_sum += tensor.gram_area(v1, v2)

        return float(dual_area_sum)


class DeficitComputer:
    """Regge deficit angle at a hinge (edge)."""

    @staticmethod
    def _dihedral_at_hinge(
        cell: IncidentCell,
        hinge: Tuple[int, int],
        coords_provider,
        tensor: MetricTensor,
    ) -> float:
        """Dihedral angle inside *cell* along the hinge edge *hinge*.

        Computed by projecting the two opposite-vertex vectors onto the plane
        orthogonal to the hinge direction and measuring the angle between them.
        """
        a, b = hinge
        others = [v for v in cell.simplex if v not in (a, b)]
        if len(others) != 2:
            raise GeometryError(
                f"Hinge {hinge} incompatible with simplex {cell.simplex}."
            )
        pa = coords_provider.vertex(a)
        pb = coords_provider.vertex(b)
        pc = coords_provider.vertex(others[0])
        pd = coords_provider.vertex(others[1])

        edge_vec = pb - pa
        vc = tensor.project_orthogonal(pc - pa, edge_vec)
        vd = tensor.project_orthogonal(pd - pa, edge_vec)

        nc = tensor.norm(vc)
        nd = tensor.norm(vd)
        if nc < 1e-12 or nd < 1e-12:
            return 0.0

        cos_theta = max(-1.0, min(1.0, tensor.inner(vc, vd) / (nc * nd)))
        return float(np.arccos(cos_theta))

    @classmethod
    def compute_hinge_deficit(
        cls,
        hinge: Tuple[int, int],
        incident_cells: Tuple[IncidentCell, ...],
        coords_provider,
        tensor: MetricTensor,
    ) -> float:
        """Regge deficit δ = 2π − Σ θ_i at the given hinge."""
        return float(
            2.0 * np.pi
            - sum(
                cls._dihedral_at_hinge(cell, hinge, coords_provider, tensor)
                for cell in incident_cells
            )
        )


# ─────────────────────────────────────────────────────────────────────────────
# 8. Engine — low-level orchestrator
# ─────────────────────────────────────────────────────────────────────────────

class GeometryEngine:
    """Low-level Regge geometry engine.

    Parameters
    ----------
    certificate    : object exposing ``tetrahedra`` (iterable of 4-tuples).
    coords_provider: object exposing ``vertex(i) → np.ndarray``.
    metric         : MetricType (optional, kept for backward compatibility).
    tensor         : MetricTensor (defaults to Euclidean of appropriate dim).
    """

    def __init__(
        self,
        certificate,
        coords_provider,
        metric: MetricType = MetricType.EUCLIDEAN,
        tensor: Optional[MetricTensor] = None,
    ) -> None:
        self.certificate = certificate
        self.coords = coords_provider
        self.metric = metric
        self._dual = BarycentricDualCell()
        # Infer ambient dimension from the first vertex if possible
        if tensor is not None:
            self.tensor = tensor
        else:
            try:
                first_tet = next(iter(certificate.tetrahedra))
                first_v = coords_provider.vertex(first_tet[0])
                dim = len(first_v)
            except Exception:
                dim = 3
            self.tensor = (
                MetricTensor.minkowski(dim)
                if metric == MetricType.LORENTZIAN
                else MetricTensor.euclidean(dim)
            )

    def run(self) -> ReggeGeometryReport:
        """Run the full geometry pipeline and return a :class:`ReggeGeometryReport`."""
        from itertools import combinations

        # ── volume + statistics pass ──────────────────────────────────────────
        total_volume = 0.0
        num_simplices = 0
        num_ignored = 0
        edge_set: set = set()
        edge_lengths: List[float] = []
        gram_volumes: List[float] = []

        for tet in self.certificate.tetrahedra:
            tet = tuple(tet)
            num_simplices += 1
            coords = tuple(self.coords.vertex(v) for v in tet)
            vol = SimplexVolume.compute_volume(
                coords, self.tensor, self.metric, strict=False
            )
            if not vol.is_valid:
                num_ignored += 1
                continue
            total_volume += vol.volume
            gram_volumes.append(vol.volume)
            for a, b in combinations(sorted(tet), 2):
                if (a, b) not in edge_set:
                    edge_set.add((a, b))
                    edge_lengths.append(
                        self.tensor.norm(
                            self.coords.vertex(b) - self.coords.vertex(a)
                        )
                    )

        stats = GeometryStatistics(
            total_volume=total_volume,
            edge_lengths=tuple(edge_lengths),
            gram_volumes=tuple(gram_volumes),
            num_simplices=num_simplices,
            num_degenerate_simplices_ignored=num_ignored,
        )

        # ── hinge extraction pass ─────────────────────────────────────────────
        hinge_to_cells: Dict[Tuple[int, int], List[IncidentCell]] = {}
        for tet in self.certificate.tetrahedra:
            tet_t = tuple(tet)
            for a, b in combinations(sorted(tet_t), 2):
                faces = tuple(
                    tuple(sorted(f))
                    for f in combinations(tet_t, 3)
                    if a in f and b in f
                )
                opposite = tuple(v for v in tet_t if v not in (a, b))
                cell = IncidentCell(tet_t, OrientationSign.POSITIVE, faces, opposite)
                hinge_to_cells.setdefault((a, b), []).append(cell)

        # ── curvature pass ────────────────────────────────────────────────────
        reports: List[CurvatureReport] = []
        for h, cells in sorted(hinge_to_cells.items()):
            cells_t = tuple(cells)
            topo = HingeTopologyAnalyzer.classify_link(cells_t)
            measure = self._dual.compute_measure(h, cells_t, self.coords, self.tensor)
            deficit = DeficitComputer.compute_hinge_deficit(
                h, cells_t, self.coords, self.tensor
            )
            scalar_curv = deficit / measure if measure > 0 else 0.0
            reports.append(
                CurvatureReport(h, cells_t, deficit, measure, measure, scalar_curv, topo)
            )

        reports_t = tuple(reports)
        num_h = len(reports_t)
        tot_def = sum(r.deficit_angle for r in reports_t)
        bulk_sum = sum(
            r.hinge_measure * r.deficit_angle
            for r in reports_t
            if r.topology == HingeTopology.INTERIOR
        )
        bdy_sum = sum(
            r.hinge_measure * r.deficit_angle
            for r in reports_t
            if r.topology == HingeTopology.BOUNDARY
        )
        total_sum = sum(r.hinge_measure * r.deficit_angle for r in reports_t)

        return ReggeGeometryReport(
            status=GeometryStatus.COMPUTED,
            total_volume=stats.total_volume,
            bulk_regge_sum=bulk_sum,
            boundary_hinge_sum=bdy_sum,
            total_geometric_sum=total_sum,
            average_scalar_curvature=tot_def / num_h if num_h > 0 else 0.0,
            maximum_deficit=max((r.deficit_angle for r in reports_t), default=0.0),
            minimum_deficit=min((r.deficit_angle for r in reports_t), default=0.0),
            total_deficit=tot_def,
            max_edge_length=max(stats.edge_lengths, default=0.0),
            min_edge_length=min(stats.edge_lengths, default=0.0),
            num_hinges_analyzed=num_h,
            num_simplices=stats.num_simplices,
            num_degenerate_simplices_ignored=stats.num_degenerate_simplices_ignored,
            min_gram_volume=min(stats.gram_volumes, default=0.0),
            max_gram_volume=max(stats.gram_volumes, default=0.0),
            curvature_reports=reports_t,
        )


# ─────────────────────────────────────────────────────────────────────────────
# 9. Universe — high-level UX entry point  (v2.6.1)
# ─────────────────────────────────────────────────────────────────────────────

class Universe:
    """High-level interface for a discrete simplicial universe.

    A *Universe* is a piecewise-linear 3-manifold (or manifold-with-boundary)
    described by a triangulation together with vertex coordinates.  It wraps
    :class:`GeometryEngine` with a convenient, stable API oriented towards
    Regge calculus and discrete quantum gravity.

    Parameters
    ----------
    tetrahedra : iterable of 4-tuples of vertex indices.
    vertices   : dict mapping vertex index → coordinate array.
    tensor     : optional :class:`MetricTensor` (default: Euclidean).

    Example
    -------
    >>> import numpy as np
    >>> verts = {0:[0,0,0], 1:[1,0,0], 2:[0.5,0.866,0], 3:[0.5,0.289,0.816]}
    >>> u = Universe([(0,1,2,3)], verts)
    >>> print(u.summary())
    """

    def __init__(self, tetrahedra, vertices: Dict, tensor: Optional[MetricTensor] = None) -> None:
        self._tetrahedra: Tuple[Tuple[int, ...], ...] = tuple(tuple(t) for t in tetrahedra)
        self._vertices: Dict[int, np.ndarray] = {
            k: np.asarray(v, dtype=float) for k, v in vertices.items()
        }
        if tensor is not None:
            self._tensor = tensor
        else:
            dim = len(next(iter(self._vertices.values()))) if self._vertices else 3
            self._tensor = MetricTensor.euclidean(dim)
        self._report: Optional[ReggeGeometryReport] = None

    # ── internals ─────────────────────────────────────────────────────────────

    class _Certificate:
        def __init__(self, tets: Tuple[Tuple[int, ...], ...]) -> None:
            self.tetrahedra = tets

    class _Coords:
        def __init__(self, v: Dict[int, np.ndarray]) -> None:
            self._v = v

        def vertex(self, i: int) -> np.ndarray:
            return self._v[i]

    # ── public API ────────────────────────────────────────────────────────────

    def run(self) -> ReggeGeometryReport:
        """Compute (or return cached) the Regge geometry of this universe."""
        if self._report is None:
            engine = GeometryEngine(
                self._Certificate(self._tetrahedra),
                self._Coords(self._vertices),
                tensor=self._tensor,
            )
            self._report = engine.run()
        return self._report

    @property
    def report(self) -> ReggeGeometryReport:
        """Cached :class:`ReggeGeometryReport` (calls :meth:`run` if needed)."""
        return self.run()

    # ── cosmological observables ──────────────────────────────────────────────

    def is_closed(self) -> bool:
        """True if every hinge is topologically interior (closed manifold)."""
        return all(
            r.topology == HingeTopology.INTERIOR for r in self.run().curvature_reports
        )

    def total_volume(self) -> float:
        """Total simplicial volume of the universe."""
        return self.run().total_volume

    def regge_action(self) -> float:
        """Total Regge–Einstein action S = Σ_h L_h δ_h."""
        return self.run().total_geometric_sum

    def bulk_action(self) -> float:
        """Regge action restricted to interior hinges."""
        return self.run().bulk_regge_sum

    def boundary_action(self) -> float:
        """Gibbons–Hawking-like contribution from boundary hinges."""
        return self.run().boundary_hinge_sum

    def topology_distribution(self) -> Dict[str, int]:
        """Count hinges by :class:`HingeTopology` type."""
        dist: Dict[str, int] = {t.name: 0 for t in HingeTopology}
        for r in self.run().curvature_reports:
            dist[r.topology.name] += 1
        return dist

    def curvature_at(self, edge: Tuple[int, int]) -> Optional[float]:
        """Deficit angle at the given edge, or *None* if not found."""
        r = self.run().get_report(edge)
        return r.deficit_angle if r is not None else None

    def summary(self) -> str:
        """Human-readable summary of the discrete universe."""
        rpt = self.run()
        dist = self.topology_distribution()
        closed = self.is_closed()
        return (
            f"=== PSGE II Universe (v{__version__}) ===\n"
            f" Topology      : {'Closed — no boundary' if closed else 'Open — has boundary'}\n"
            f" Simplices     : {rpt.num_simplices}"
            f" ({rpt.num_degenerate_simplices_ignored} degenerate ignored)\n"
            f" Hinges        : {rpt.num_hinges_analyzed}"
            f" (Interior {dist.get('INTERIOR', 0)},"
            f" Boundary {dist.get('BOUNDARY', 0)},"
            f" Non-manifold {dist.get('NON_MANIFOLD', 0)})\n"
            f" Total Volume  : {rpt.total_volume:.6f}\n"
            f" Regge Action  : {rpt.total_geometric_sum:.6f}"
            f"  [bulk {rpt.bulk_regge_sum:.6f}"
            f" | boundary {rpt.boundary_hinge_sum:.6f}]\n"
            f" Curvature     : avg {rpt.average_scalar_curvature:.4f} rad/unit,"
            f" range [{rpt.minimum_deficit:.4f}, {rpt.maximum_deficit:.4f}] rad\n"
            f" Edge lengths  : min {rpt.min_edge_length:.4f},"
            f" max {rpt.max_edge_length:.4f}\n"
            f"{'=' * 42}"
        )
