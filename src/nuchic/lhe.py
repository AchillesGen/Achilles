""" Implement a LHE writer. """

import xml.etree.ElementTree as ElementTree
import xml.dom.minidom as minidom

# TODO: Improve documentation
# TODO: Implement additional version specific features


def cast(lst, dtype):
    """ Convert a list to a list of ints. """
    if isinstance(lst, list):
        return list(map(dtype, lst))
    return dtype(lst)


class LHEFExcept(Exception):
    """ Exception class for reading LHE tags. """
    def __init__(self, element, tag, expected):
        message = 'LHE: Cannot construct {} from: {} expected {}'.format(
            element, tag, expected)
        super(LHEFExcept, self).__init__(message)
        self.message = message


class LHAProcess:
    """ Store information about the processes in the LHE file. """
    def __init__(self, xsec, xerr, xmax, processid):
        self.xsec = cast(xsec, float)
        self.xerr = cast(xerr, float)
        self.xmax = cast(xmax, float)
        self.proc_id = cast(processid, int)

    def __str__(self):
        return ' {:>14} {:>14} {:>14} {:>6}\n'.format(
            self.xsec, self.xerr, self.xmax, self.proc_id)

    def __repr__(self):
        return 'LHAProcess({}, {}, {}, {})'.format(
            self.xsec, self.xerr, self.xmax, self.proc_id)

    def __eq__(self, other):
        return (self.xsec == other.xsec and
                self.xerr == other.xerr and
                self.xmax == other.xmax and
                self.proc_id == other.proc_id)


class LHAParticle:
    """ Store information about the particle in the LHE file. """
    def __init__(self, pid, status, momentum, **kwargs):
        self.pid = cast(pid, int)
        self.status = cast(status, int)
        self.mom = cast(momentum, float)
        self.col = cast(kwargs['col'], int) if 'col' in kwargs else [0, 0]
        self.mothers = cast(kwargs['mothers'], int) \
            if 'mothers' in kwargs else [0, 0]
        self.vtim = cast(kwargs['vtim'], float) if 'vtim' in kwargs else 0
        self.spin = cast(kwargs['spin'], int) if 'spin' in kwargs else 9

    def __str__(self):
        return (' {:>8} {:>2} {:>4} {:>4} {:>4} {:>4} {:>18} {:>18} {:>18} '
                '{:>18} {:>18} {:>1} {:>1}\n'.format(self.pid, self.status,
                                                     *self.mothers, *self.col,
                                                     *self.mom, self.vtim,
                                                     self.spin))

    def __repr__(self):
        return 'LHAParticle({}, {}, {}, {}, {}, {}, {})'.format(
            self.pid, self.status, self.mothers, self.col,
            self.mom, self.vtim, self.spin)

    def __eq__(self, other):
        return (self.pid == other.pid and
                self.status == other.status and
                self.mom == other.mom and
                self.col == other.col and
                self.mothers == other.mothers and
                self.vtim == other.vtim and
                self.spin == other.spin)


class TagFactory:
    """ Factory class to create the right tag. """
    def __init__(self):
        self._tags = {}

    def register_tag(self, key, tag):
        """ Register a tag to the factory. """
        self._tags[key] = tag

    def create(self, key):
        """ Create tag from a key. """
        tag = self._tags.get(key)
        if not tag:
            return LHATag(key)
        return tag()

    def __call__(self, key):
        return self.create(key)


LHAFACTORY = TagFactory()


class LHATag:
    """ Base class for LHATags. """
    def __init__(self, tagname):
        self.tagname = tagname
        self.attrib = {}
        self.children_names = []
        self.children = {}
        self.params = {}
        self.text = ''

    @property
    def name(self):
        """ Return name of the tag. """
        name = ''
        if 'name' in self.attrib:
            name = self.attrib['name']
        elif 'id' in self.attrib:
            name = self.attrib['id']
        else:
            name = self.tagname
        return name

    def from_xml(self, element):
        """ Read the element from a file. """
        # Ensure it has the proper tag
        if element.tag != self.tagname:
            raise LHEFExcept(self.__class__.__name__,
                             element.tag,
                             self.tagname)

        # Find all keys
        for key, value in element.attrib.items():
            self.attrib[key] = value

        # Find all children
        for child_name in self.children_names:
            for child in element.findall(child_name):
                child_element = LHAFACTORY(child_name)
                child_element.from_xml(child)
                name = '{}.{}'.format(child_element.tagname,
                                      child_element.name)
                self.children[name] = child_element

        self.text = element.text if element.text is not None else ''

    def to_xml(self, root):
        """ Write the element to a file. """
        xml = ElementTree.SubElement(root, self.tagname)
        for key, value in self.attrib.items():
            xml.attrib[key] = str(value)
        for key, value in self.children.items():
            value.to_xml(xml)
        xml.text = str(self.text)

    def __eq__(self, other):
        """ Method to compare two tags are equal. """
        return \
            (self.name == other.name and
             self.attrib == other.attrib and
             len(self.children_names) == len(other.children_names) and
             len(self.children) == len(other.children) and
             self.params == other.params and
             self.text == other.text)


class LHAweights(LHATag):
    """ Store weights classes. """
    def __init__(self):
        super().__init__('weights')

    @property
    def weights(self):
        """ Access to the weights stored in text. """
        return self.text.split(' ')

    def add_weights(self, weights):
        """ Add weights to the tag. """
        self.text = ' '.join([str(tmp) for tmp in weights])


LHAFACTORY.register_tag('weights', LHAweights)


class LHAscales(LHATag):
    """ LHA tag to store scales of the calculation. """
    def __init__(self):
        super().__init__('scales')

    @property
    def mur(self):
        """ Get renormalization scale. """
        return float(self.attrib['mur'])

    @mur.setter
    def mur(self, mur):
        self.attrib['mur'] = mur

    @property
    def muf(self):
        """ Get factorization scale. """
        return float(self.attrib['muf'])

    @muf.setter
    def muf(self, muf):
        self.attrib['muf'] = muf

    @property
    def mups(self):
        """ Get parton shower scale. """
        return float(self.attrib['mups'])

    @mups.setter
    def mups(self, mups):
        self.attrib['mups'] = mups


LHAFACTORY.register_tag('scales', LHAscales)


class LHAgenerator(LHATag):
    """ LHA tag storing information about the generator. """
    def __init__(self):
        super().__init__('generator')

    def create(self, name, version, comments='', attrib=None):
        """ Create needed info for a generator. """
        self.attrib['name'] = name
        self.attrib['version'] = version
        self.text = comments
        if attrib is not None:
            self.attrib = attrib


LHAFACTORY.register_tag('generator', LHAgenerator)


class LHAwgt(LHATag):
    """ LHA tag storing information about the weight as text. """
    def __init__(self):
        super().__init__('wgt')

    @property
    def identity(self):
        """ Access to the id of the wgt. """
        return self.attrib['id']

    @identity.setter
    def identity(self, identity):
        self.attrib['id'] = str(identity)

    @property
    def weight(self):
        """ Access to the value of the wgt. """
        return float(self.text)

    @weight.setter
    def weight(self, weight):
        self.text = str(weight)


LHAFACTORY.register_tag('wgt', LHAwgt)


class LHAweight(LHATag):
    """ LHA tag storing information about the weight. """
    def __init__(self):
        super().__init__('weight')

    @property
    def identity(self):
        """ Access to the id of the wgt. """
        return self.attrib['id']

    @identity.setter
    def identity(self, identity):
        self.attrib['id'] = str(identity)

    @property
    def weights(self):
        """ Access to the weight dictionary. """
        weights = {}
        for weight in self.text.split('  '):
            key, value = weight.split(' = ')
            weights[key] = float(value)
        return weights

    def add_weight(self, name, value):
        """ Add a weight. """
        self.text += ' {} = {} '.format(name, value)


LHAFACTORY.register_tag('weight', LHAweight)


class LHAweightgroup(LHATag):
    """ LHA tag storing all the different weight information. """
    def __init__(self):
        super().__init__('weightgroup')
        self.children_names = ['weight']

    def add_weight(self, weights):
        """ Add a single weight's information to the group. """
        weight = LHAweight()
        weight.identity = len(self.children)+1
        for key, value in weights.items():
            if key == 'id':
                weight.identity = str(value)
            else:
                weight.add_weight(key, value)
        self.children['weight.' + weight.identity] = weight

    def add_weights(self, weights):
        """ Add multiple weights' information to the group. """
        for weight in weights:
            self.add_weight(weight)

    @property
    def name(self):
        """ Get name of group from type information. """
        if 'type' in self.attrib:
            return self.attrib['type']
        return super().name


LHAFACTORY.register_tag('weightgroup', LHAweightgroup)


class LHArwgt(LHATag):
    """ LHA Tag for reweighting wgts. """
    def __init__(self):
        super().__init__('rwgt')
        self.children_names = ['wgt']

    def add_weight(self, identity, wgt):
        """ Add a reweighting weight. """
        tmp = LHAwgt()
        tmp.identity = identity
        tmp.weight = wgt
        self.children['wgt.' + str(identity)] = tmp


LHAFACTORY.register_tag('rwgt', LHArwgt)


class LHAinitrwgt(LHATag):
    """ LHA Tag to initialize reweighting block. """
    def __init__(self):
        super().__init__('initrwgt')
        self.children_names = ['weightgroup', 'weight']

    def add_weightgroup(self, weightgroup):
        """ Add a weight group to reweighting header. """
        self.children['weightgroup.' + weightgroup.name] = weightgroup

    def add_weight(self, weight):
        """ Add a weight to reweighting header. """
        self.children['weight.' + weight.name] = weight


LHAFACTORY.register_tag('initrwgt', LHAinitrwgt)


# TODO: Pretty up this class
class LHAinit(LHATag):
    """ LHA Tag for initializing the process block. """
    def __init__(self):
        super().__init__('init')
        self.children_names = ['generator']
        self.params['beams'] = []
        self.params['energy'] = []
        self.params['pdfgroup'] = []
        self.params['pdfset'] = []
        self.params['wgtid'] = -4
        self.params['processes'] = []

    @property
    def beams(self):
        """ The initial state beam constituents. """
        return self.params['beams']

    @beams.setter
    def beams(self, beams):
        self.params['beams'] = beams

    @property
    def energy(self):
        """ The energy of each beam. """
        return self.params['energy']

    @energy.setter
    def energy(self, energy):
        self.params['energy'] = energy

    @property
    def pdfgroup(self):
        """ The pdf group for each of the beams. """
        return self.params['pdfgroup']

    @pdfgroup.setter
    def pdfgroup(self, pdfgroup):
        self.params['pdfgroup'] = pdfgroup

    @property
    def pdfset(self):
        """ The pdf set for each of the beams. """
        return self.params['pdfset']

    @pdfset.setter
    def pdfset(self, pdfset):
        self.params['pdfset'] = pdfset

    @property
    def wgtid(self):
        """ The definition of what the weight corresponds to. """
        return self.params['wgtid']

    @wgtid.setter
    def wgtid(self, wgtid):
        self.params['wgtid'] = wgtid

    @property
    def processes(self):
        """ A list of all the processes contained within the LHE file. """
        return self.params['processes']

    def add_process(self, process):
        """ Add a process to the list of processes. """
        self.params['processes'].append(process)

    def from_xml(self, element):
        super().from_xml(element)
        # Split text by lines
        lines = [s.rstrip() for s in self.text.splitlines() if s.strip()]

        # Load basic details into params
        init = lines[0].split()
        self.params['beams'] = cast(init[:2], int)
        self.params['energy'] = cast(init[2:4], int)
        self.params['pdfgroup'] = cast(init[4:6], int)
        self.params['pdfset'] = cast(init[6:8], int)
        self.params['wgtid'] = cast(init[8], int)
        nprocesses = cast(init[9], int)

        # Load information for each process
        for line in lines[1:1+nprocesses]:
            tmp = line.split()
            self.params['processes'].append(
                LHAProcess(*tmp))
        self.text = '\n'.join(lines[1+nprocesses:]).lstrip()

    def to_xml(self, root):
        text = '\n {:>8} {:>8} '.format(*self.beams)
        text += '{:>14} {:>14} '.format(*self.energy)
        text += '{:>4} {:>4} '.format(*self.pdfgroup)
        text += '{:>4} {:>4} '.format(*self.pdfset)
        text += '{:>4} {:>4}\n'.format(self.wgtid,
                                       len(self.params['processes']))
        for process in self.params['processes']:
            text += str(process)
        tmp = self.text
        self.text = (text + '\t' + self.text).rstrip()
        super().to_xml(root)
        self.text = tmp


LHAFACTORY.register_tag('init', LHAinit)


# TODO: Pretty up this class
class LHAevent(LHATag):
    """ LHA Tag to build an event block. """
    def __init__(self):
        super().__init__('event')
        self.children_names = ['rwgt', 'weights', 'scales']
        self.params['particles'] = []
        self.params['nprocess'] = -1
        self.params['wgt'] = 0
        self.params['scale'] = 0
        self.params['aqed'] = 0
        self.params['aqcd'] = 0

    @property
    def particles(self):
        """ Get a list of all particles. """
        return self.params['particles']

    def add_particle(self, particle):
        """ Add a particle to the event. """
        self.params['particles'].append(particle)

    @property
    def nprocess(self):
        """ The event's process id. """
        return self.params['nprocess']

    @nprocess.setter
    def nprocess(self, nprocess):
        self.params['nprocess'] = nprocess

    @property
    def wgt(self):
        """ The event's weight. """
        return self.params['wgt']

    @wgt.setter
    def wgt(self, wgt):
        self.params['wgt'] = wgt

    @property
    def scale(self):
        """ The hard process scale of the event. """
        return self.params['scale']

    @scale.setter
    def scale(self, scale):
        self.params['scale'] = scale

    @property
    def aqed(self):
        """ The value of alpha_qed at the hard scale of the event. """
        return self.params['aqed']

    @aqed.setter
    def aqed(self, aqed):
        self.params['aqed'] = aqed

    @property
    def aqcd(self):
        """ The value of alpha_qed at the hard scale of the event. """
        return self.params['aqcd']

    @aqcd.setter
    def aqcd(self, aqcd):
        self.params['aqcd'] = aqcd

    def from_xml(self, element):
        super().from_xml(element)
        lines = self.text.split('\n')[1:]
        header = lines[0].split()
        nparts = cast(header[0], int)
        self.params['nprocess'] = cast(header[1], int)
        self.params['wgt'] = cast(header[2], float)
        self.params['scale'] = cast(header[3], float)
        self.params['aqed'] = cast(header[4], float)
        self.params['aqcd'] = cast(header[5], float)
        particles = lines[1:1+nparts]

        for part in particles:
            part_info = part.split()
            # Information is given in the following format
            # - Particle id
            # - Particle status
            # - Particle mothers
            # - Particle color and anti-color
            # - Particle 4-momentum and mass
            # - Particle lifetime
            # - Particle spin mode
            # See LHAParticle for additional details
            self.params['particles'].append(
                LHAParticle(part_info[0], part_info[1],
                            [part_info[6], part_info[7],
                             part_info[8], part_info[9],
                             part_info[10]],
                            col=[part_info[4], part_info[5]],
                            mothers=[part_info[2], part_info[3]],
                            vtim=part_info[11], spin=part_info[12]))
        self.text = '\n'.join(lines[1+nparts:]).rstrip()

    def to_xml(self, root):
        text = '\n {:>4} {:>6} {:>14} {:>14} {:>14} {:>14}\n'.format(
            len(self.params['particles']), self.params['nprocess'],
            self.params['wgt'], self.params['scale'], self.params['aqed'],
            self.params['aqcd'])
        for part in self.params['particles']:
            text += str(part)
        tmp = self.text
        self.text = text + self.text + '    '
        super().to_xml(root)
        self.text = tmp


LHAFACTORY.register_tag('event', LHAevent)


class LHEF:
    """ Class to read and write LHEF files. """
    def __init__(self, filename):
        self.filename = filename
        self.events = []
        self.root = None
        self.version = None
        self.header = None
        self.init = None

    def read(self):
        """ Read an LHEF file into memory. """
        tree = ElementTree.parse(self.filename)
        self.root = tree.getroot()
        if self.root.tag != 'LesHouchesEvents':
            raise ValueError('Invalid File Format: Incorrect header, '
                             'found {}'.format(self.root.tag))
        try:
            self.version = self.root.attrib['version']
        except KeyError:
            raise KeyError('Invalid File Format: No version number found.')

        self.header = self.root.find('header')
        self.header.text = self.header.text.strip()
        for child in self.header:
            if child.tag == 'initrwgt':
                initrwgt = LHAinitrwgt()
                initrwgt.from_xml(child)

        self.init = LHAinit()
        self.init.from_xml(self.root.find('init'))

        self.events = []
        for event_xml in self.root.iter('event'):
            event = LHAevent()
            event.from_xml(event_xml)
            self.events.append(event)

    def set_header(self, version, initrwgt=None, header_comments=None):
        """ Set the header block for the LHE file. """
        self.version = str(version)
        self.root = ElementTree.Element('LesHouchesEvents')
        self.root.attrib['version'] = str(version)
        self.header = ElementTree.SubElement(self.root, 'header')
        if header_comments is not None:
            self.header.text = header_comments
        else:
            self.header.text = ''
        if initrwgt is not None:
            initrwgt.to_xml(self.header)

    def set_init_block(self, init):
        """ Set the initialization block for the LHE file. """
        self.init = init
        init.to_xml(self.root)

    def add_event(self, event):
        """ Add an event to the xml tree. """
        event.to_xml(self.root)
        self.events.append(event)

    def write(self, outfile=None):
        """ Write out the LHE file. """
        xmlstr = minidom.parseString(
            ElementTree.tostring(self.root)).toprettyxml(indent='    ')

        if outfile is None:
            with open(self.filename, 'w') as out:
                out.write(xmlstr[23:])
        else:
            with open(outfile, 'w') as out:
                out.write(xmlstr[23:])


if __name__ == '__main__':
    def main():
        """ Test function. """
        test = LHEF('test_write.lhe')
        initrwgt = LHAinitrwgt()
        weightgroup = LHAweightgroup()
        weightgroup.attrib['type'] = 'scale_variation'
        weightgroup.attrib['combine'] = 'envelope'
        weights = [{'mur': 0.1E1, 'muf': 0.1E1},
                   {'mur': 0.1E1, 'muf': 0.2E1},
                   {'mur': 0.2E1, 'muf': 0.1E1},
                   {'mur': 0.2E1, 'muf': 0.2E1},
                   {'mur': 0.1E1, 'muf': 0.5E0},
                   {'mur': 0.5E0, 'muf': 0.1E1},
                   {'mur': 0.5E0, 'muf': 0.5E0}]
        weightgroup.add_weights(weights)

        initrwgt.children[weightgroup.name] = weightgroup
        test.set_header(3.0, initrwgt)

        init = LHAinit()
        init.beams = [2212, 2212]
        init.energy = [4000, 4000]
        init.pdfgroup = [-1, -1]
        init.pdfset = [21100, 21100]
        init.wgtid = -4
        process = LHAProcess(0.501090862e2, 0.89185414e-1, 0.50109093e2, 66)
        init.processes.append(process)
        gen = LHAgenerator()
        gen.create('SomeGen1', '1.2.3', 'Some additional comments')
        init.children['generator.SomeGen1'] = gen

        test.set_init_block(init)

        event = LHAevent()
        event.attrib['npLO'] = ' -1 '
        event.attrib['npNLO'] = ' 1 '
        event.nprocess = 66
        event.wgt = 0.501
        event.scale = 0.14137
        event.aqed = 0.755e-2
        event.aqcd = 0.1211

        part = LHAParticle(5, -1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
        event.particles.append(part)

        part = LHAParticle(2, -1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
        event.particles.append(part)

        part = LHAParticle(24, 1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
        event.particles.append(part)

        part = LHAParticle(2, 1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
        event.particles.append(part)

        test.add_event(event)
        test.write()

    main()
