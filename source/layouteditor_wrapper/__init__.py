"""
Create a new design::

    /some/directory/somewhere> python -i -m layouteditor_wrapper
    Variable 'l' is <class 'layouteditor_wrapper.wrapper.Layout'>
    Variable 'd' is <class 'layouteditor_wrapper.wrapper.Drawing'>
    >>> d.cells
    OrderedDict([('noname', <layouteditor_wrapper.wrapper.Cell object at 0x000001CDFF756B38>)])

Open an existing design::

    /some/directory/somewhere> python -i -m layouteditor_wrapper existing_design.gds
    Variable 'l' is <class 'layouteditor_wrapper.wrapper.Layout'>
    Variable 'd' is <class 'layouteditor_wrapper.wrapper.Drawing'>
    >>> d.cells.keys()
    odict_keys([<cell names in the existing design>])

"""
