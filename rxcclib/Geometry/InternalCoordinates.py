
class InternalCoordinates(object):
    def __init__(self, fileobj, atoms):
        self.file = fileobj
        assert isinstance(atoms, list), 'give atoms as list'
        self.atomset = atoms
        if len(atoms) == 2:
            self.type = 'stretching'
        elif len(atoms) == 3:
            self.type = 'bending'
        elif len(atoms) == 4:
            self.type = 'torsion'
        else:
            raise


