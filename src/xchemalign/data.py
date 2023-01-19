import gemmi
from pydantic import BaseModel, validator

Structure = gemmi.Structure


class LigandID(BaseModel):
    dtag: str
    id: int

    def __eq__(self, other) -> bool:
        if self.dtag == other.dtag:
            if self.id == other.id:
                return True
        return False

    def __hash__(self):
        return hash((self.dtag, self.id))


class SymOp(BaseModel):
    operation: str
    image: tuple[int, int, int]


class AtomID(BaseModel):
    chain: str
    residue: int
    atom: str

    def __eq__(self, other) -> bool:
        if self.atom == other.atom:
            if self.residue == other.residue:
                if self.chain == other.chain:
                    return True

        return False

    def __hash__(self):
        return hash((self.chain, self.residue, self.atom))


class Atom(BaseModel):
    element: str
    atom_id: AtomID
    x: float
    y: float
    z: float


class AlignableSite(BaseModel):
    id: int
    name: str
    ligand_ids: list[LigandID]


class CanonicalSite(BaseModel):
    id: int
    name: str
    refpdb: str
    atoms: dict[AtomID, Atom]
    literatureref: str
    members: list[LigandID]


class XtalForm(BaseModel):
    id: int
    space_group: int
    unit_cell: tuple[float, float, float, float, float, float]


class XtalFormSite(BaseModel):
    id: int
    canon_site_id: int
    xtal_form_id: int
    code: str
    refpdb: str
    atoms: dict[AtomID, Atom]
    artefact_atoms: dict[AtomID, Atom]
    members: list[LigandID]


class SiteObservation(BaseModel):
    id: int
    ligand_id: LigandID
    xtal_form_site_id: int
    fragalysis_label: str
    description: str
    code: str
    dataset: str
    # compound: int


class Dataset(BaseModel):
    dtag: str
    pdb: str
    # mtz: str


class SystemData(BaseModel):
    dataset: list[Dataset]


class LigandNeighbourhood(BaseModel):
    atoms: dict[AtomID, Atom]
    artefact_atoms: dict[AtomID, Atom]


class DatasetID:
    dtag: str


class SystemSites(BaseModel):
    canonical_site: dict[int, CanonicalSite]
    xtal_form_site: dict[int, XtalFormSite]
    site_observation: dict[LigandID, SiteObservation]

    @validator("canonical_site")
    def check_canonical_site_ids(cls, v: dict[LigandID, SiteObservation]):
        if not v:
            return
        for site_id, site in v.items():
            assert site_id == site.id

    @validator("canonical_site")
    def check_canonical_site_ids_sequential(
        cls,
        v: dict[LigandID, SiteObservation],
    ):
        if not v:
            return
        num_sites: int = len(v)
        site_nums = [site.id for site in v.values()]
        for site_num in range(num_sites):
            assert site_num in site_nums

    @validator("xtal_form_site")
    def check_xtal_form_site_ids(cls, v: dict[int, XtalFormSite]):
        if not v:
            return
        for site_id, site in v.items():
            assert site_id == site.id

    @validator("xtal_form_site")
    def check_xtal_form_site_ids_sequential(
        cls,
        v: dict[int, XtalFormSite],
    ):
        if not v:
            return
        num_sites: int = len(v)
        site_nums = [site.id for site in v.values()]
        for site_num in range(num_sites):
            assert site_num in site_nums

    @validator("site_observation")
    def check_site_observation_ids(cls, v: dict[LigandID, SiteObservation]):
        if not v:
            return
        for site_id, site in v.items():
            assert site_id == site.ligand_id

    @validator("site_observation")
    def check_site_observation_ids_sequential(
        cls,
        v: dict[LigandID, SiteObservation],
    ):
        if not v:
            return
        num_sites: int = len(v)
        site_nums = [site.id for site in v.values()]
        for site_num in range(num_sites):
            assert site_num in site_nums
