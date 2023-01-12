from pydantic import BaseModel


class LigandID(BaseModel):
    dtag: str
    id: int


class SymOp(BaseModel):
    operation: str
    image: tuple[int, int, int]


class AtomID(BaseModel):
    chain: str
    residue: int
    atom: str


class CanonicalSite(BaseModel):
    id: int
    name: str
    refpdb: str
    atoms: list[tuple[AtomID, SymOp]]
    literatureref: str
    members: list[LigandID]


class XtalForm(BaseModel):
    id: int
    space_group: int
    unit_cell: tuple[float, float, float, float, float, float]


class XtalFormSite(BaseModel):
    canon_site_id: int
    xtal_form_id: int
    code: str
    refpdb: str
    atoms: list[tuple[AtomID, SymOp]]
    artefact_atoms: list[tuple[AtomID, SymOp]]
    members: list[LigandID]


class SiteObservation(BaseModel):
    id: int
    xtal_form_site_id: int
    fragalysis_label: str
    description: str
    code: str
    dataset: int
    compound: int


class Dataset(BaseModel):
    dtag: str
    pdb: str
    # mtz: str


class SystemData(BaseModel):
    dataset: list[Dataset]


class LigandNeighbourhood(BaseModel):
    atoms: list[tuple[AtomID, SymOp]]
    artefact_atoms: list[tuple[AtomID, SymOp]]


class DatasetID:
    dtag: str


class SystemSites(BaseModel):
    canonical_site: dict[int, CanonicalSite]
    xtal_form_site: dict[int, XtalFormSite]
    site_observation: dict[LigandID, SiteObservation]
