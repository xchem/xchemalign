import gemmi
from pydantic import BaseModel, validator

Structure = gemmi.Structure


class DatasetID(BaseModel):
    dtag: str


class LigandID(BaseModel):
    dtag: str
    chain: str
    id: int

    def __eq__(self, other) -> bool:
        if self.dtag == other.dtag:
            if self.id == other.id:
                if self.chain == other.chain:
                    return True
        return False

    def __hash__(self):
        return hash((self.dtag, self.chain, self.id))


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


class Transform(BaseModel):
    vec: list[float]
    mat: list[list[float]]


class Transforms(BaseModel):
    ligand_ids: list[tuple[LigandID, LigandID]]
    transforms: list[Transform]

    def get_transform(self, transform_id: tuple[LigandID, LigandID]):
        for ligand_id, _transform in zip(self.ligand_ids, self.transforms):
            if transform_id == ligand_id:

                transform = gemmi.Transform()
                transform.vec.fromlist(_transform.vec)
                transform.mat.fromlist(_transform.mat)

                return transform


class Atom(BaseModel):
    element: str
    atom_id: AtomID
    x: float
    y: float
    z: float
    image: Transform


class AlignableSite(BaseModel):
    id: int
    name: str
    ligand_ids: list[LigandID]
    reference: LigandID


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
    members: list[DatasetID]
    transforms: list[Transform]


class XtalForms(BaseModel):
    xtalforms: list[XtalForm]


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


class LigandBindingEvent(BaseModel):
    id: int
    chain: str
    residue: int
    xmap: str


class LigandBindingEvents(BaseModel):
    ligand_ids: list[LigandID]
    ligand_binding_events: list[LigandBindingEvent]

    def __getitem__(self, lid: LigandID) -> LigandBindingEvent | None:
        for _lid, lbe in zip(self.ligand_ids, self.ligand_binding_events):
            if lid == _lid:
                return lbe

        return None


class Dataset(BaseModel):
    dtag: str
    pdb: str
    ligand_binding_events: LigandBindingEvents
    # mtz: str


class SystemData(BaseModel):
    dataset_ids: list[DatasetID]
    datasets: list[Dataset]

    def get_dataset(self, did: str | DatasetID) -> Dataset:
        for _did, dataset in zip(self.dataset_ids, self.datasets):
            if _did == did:
                return dataset

        raise Exception()


class LigandNeighbourhood(BaseModel):
    atom_ids: list[AtomID]
    atoms: list[Atom]
    artefact_atom_ids: list[AtomID]
    artefact_atoms: list[Atom]


class LigandNeighbourhoods(BaseModel):
    ligand_ids: list[LigandID]
    ligand_neighbourhoods: list[LigandNeighbourhood]

    def get_neighbourhood(self, ligand_id: LigandID):
        for _ligand_id, _neighbourhood in zip(
            self.ligand_ids, self.ligand_neighbourhoods
        ):
            if _ligand_id == ligand_id:
                return _neighbourhood


class SystemSites(BaseModel):
    canonical_site: dict[int, CanonicalSite]
    xtal_form_site: dict[int, XtalFormSite]
    ligand_ids: list[LigandID]
    site_observation: dict[int, SiteObservation]

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


class Site(BaseModel):
    members: list[LigandID]


class Sites(BaseModel):
    site_ids: list[int]
    sites: list[AlignableSite]


class Block(BaseModel):
    xi: int
    yi: int
    zi: int
    dx: int
    dy: int
    dz: int
    transform: Transform
