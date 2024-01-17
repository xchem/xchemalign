XChemAlign again 
===========================

|code_ci| |docs_ci| |coverage| |pypi_version| |license|


XChemAlign is a small command line interface and library for constructing data models of fragment screening results, identifying groupings of fragments and aligning structures and ccp4 maps into a common reference frame.

XChemAlign is designed to work in the presence of crystallographic symmetry, non-crystallographic symmetry, multiple space groups and conformational hetrogeneity among the source data.

If you have seen XChemAlign fail to correctly parse data, identify sites or subsites, align structures or xmaps or otherwise have any suggestions, please raise an issue on this github!

============== ==============================================================
PyPI           ``pip install xchemalign``
Source code    https://github.com/xchem/xchemalign
Documentation  https://xchem.github.io/xchemalign
Releases       https://github.com/xchem/xchemalign/releases
============== ==============================================================

The best way to install is to create an anaconda environment and then run::
    git clone 

The prefferred and more reproducible way to use XChemAlign to create a new project is to create an options json::

    {
        "source_dir": "/path/to/source/dir", 
        "datasources": ["/path/to/datasource/1", "/path/to/datasource/2"], 
        "panddas": ["/path/to/pandda/1", "/path/to/pandda/2"]
    }

And then run the process option. This will produce aligned structures and xmaps in the source_dir::

    $ python -m ligand_neighbourhood_alignment.cli process /path/to/option/json

However, XChemAlign can be manually used as follows. This begins with initializing a project::

    $ python -m ligand_neighbourhood_alignment.cli init /path/to/project/dir

After this datasources can be added::

    $ python -m ligand_neighbourhood_alignment.cli add_datasource /path/to/project/dir /path/to/datasource

In order to align PanDDA event maps, the PanDDA directory in which these maps were modelled with PanDDA inspect must be added::

    $ python -m ligand_neighbourhood_alignment.cli add_datasource /path/to/project/dir /path/to/pandda

Multiple datasources can be added to a single project, as can multiple PanDDAs. Be cautious if the fragment was modelled in multiple of these PanDDAs, as only the xmaps from the first added PanDDA will be included.

After sources and PanDDAs have been added, they can be parsed to produced a combined dataset::

    $ python -m ligand_neighbourhood_alignment.cli parse_data_sources /path/to/project/dir

Once the combined datasource has been produced, the alignment graph can be calculated. This graph encodes the transformations necessary to align any two fragments based on their local protein environment (assuming this is possible)::

    $ python -m ligand_neighbourhood_alignment.cli build_graph /path/to/project/dir

Once the local alignment graph has been determined, subsites can be calculated based on the connected components (ligands that can be directly or indirectly aligned to one another based on their local protein environment), and then sites can be determined from those subsites which share significant numbers of residues::

    $ python -m ligand_neighbourhood_alignment.cli generate_sites_from_components /path/to/project/dir

Once the sites and subsites have been generated, all of the structures can be aligned to a common reference frame. This first aligns all structures to the reference ligand in their subsite, then all subsites in a site to their reference subsite, then all sites to the reference site:: 

    $ python -m ligand_neighbourhood_alignment.cli align_structures /path/to/project/dir

Finally, the refined 2Fo-Fc maps and any event maps from PanDDAs can be aligned to their respective aligned structures::

    $ python -m ligand_neighbourhood_alignment.cli align_xmaps /path/to/project/dir

In order to inspect all of the structures in a site, there is a convenience function to open them all at once in coot::

    $ python -m ligand_neighbourhood_alignment.cli open_site /path/to/project/dir site_number

If new datasources or PanDDAs have been added, or manual changes have been made to the sites or subsites, then these changes can be conveniently propogated to a new set of final sites, subsites, structures and xmaps using the update command::

    $ python -m ligand_neighbourhood_alignment.cli update /path/to/project/dir


.. |code_ci| image:: https://github.com/ConorFWild/ligand_neighbourhood_alignment/actions/workflows/code.yml/badge.svg?branch=main
    :target: https://github.com/ConorFWild/ligand_neighbourhood_alignment/actions/workflows/code.yml
    :alt: Code CI

.. |docs_ci| image:: https://github.com/ConorFWild/ligand_neighbourhood_alignment/actions/workflows/docs.yml/badge.svg?branch=main
    :target: https://github.com/ConorFWild/ligand_neighbourhood_alignment/actions/workflows/docs.yml
    :alt: Docs CI

.. |coverage| image:: https://codecov.io/gh/xchem/xchemalign/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/ConorFWild/ligand_neighbourhood_alignment
    :alt: Test Coverage

.. |pypi_version| image:: https://img.shields.io/pypi/v/xchemalign.svg
    :target: https://pypi.org/project/xchemalign
    :alt: Latest PyPI version

.. |license| image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://opensource.org/licenses/Apache-2.0
    :alt: Apache License

..
    Anything below this line is used when viewing README.rst and will be replaced
    when included in index.rst

See https://xchem.github.io/xchemalign for more detailed documentation.
