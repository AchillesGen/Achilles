""" Testing the LHE reader and writer. """

# TODO: Reduce amount of repeated code by adding fixtures

from xml.etree import ElementTree
import numpy as np

import nuchic.lhe as lhe


def test_init_lha_tag():
    """ Test LHATags base class. """
    tag = lhe.LHATag('test')
    assert tag.name == 'test'
    assert tag.attrib == {}
    assert tag.children_names == []
    assert tag.children == {}
    assert tag.params == {}
    assert tag.text == ''


def test_lha_tag_to_from_xml_basic():
    """ Test LHATags basic to and from xml functions. """
    tag = lhe.LHATag('test')
    tag.attrib['foo'] = 'bar'
    tag.text = 'this is some test text.'

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHATag('test')

    head = root.findall('test')
    for child in head:
        assert child.tag == 'test'
        tag2.from_xml(child)

    assert tag == tag2


def test_lha_tag_with_child():
    """ Test LHATags with child to and from xml functions. """
    tag_child = lhe.LHATag('child')
    tag_child.attrib['baz'] = 'foo'
    tag_child.text = 'this is some text from the child.'

    tag = lhe.LHATag('parent')
    tag.attrib['foo'] = 'bar'
    tag.text = 'this is some text from the parent.'
    tag.children['child.child'] = tag_child
    tag.children_names = ['child']

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHATag('parent')
    tag2.children_names = ['child']

    head = root.findall('parent')
    for parent in head:
        assert parent.tag == 'parent'
        tag2.from_xml(parent)

    assert tag == tag2

    for key, _ in tag.children.items():
        assert tag.children[key] == tag2.children[key]


def test_lha_weights():
    """ Test the LHAweights class. """
    tag = lhe.LHAweights()
    weights = np.random.random(5)
    tag.add_weights(weights)

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHAweights()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2
    assert np.allclose(weights, np.array(tag2.weights, dtype=np.float64))


def test_lha_scales():
    """ Test the LHAscales class. """
    tag = lhe.LHAscales()
    tag.mur = '1'
    tag.muf = '1'
    tag.mups = '1'

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHAscales()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2
    assert tag.mur == tag2.mur
    assert tag.muf == tag2.muf
    assert tag.mups == tag2.mups


def test_lha_generator():
    """ Test LHA Generator tag. """
    tag = lhe.LHAgenerator()
    tag.create('SomeGen1', '1.2.3', 'Some additional comments',
               attrib={'foo': 'bar'})

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHAgenerator()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2


def test_lha_wgt():
    """ Test LHA wgt tag. """
    tag = lhe.LHAwgt()
    tag.identity = 0
    tag.weight = 0.5

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHAwgt()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2
    assert tag.identity == tag2.identity
    assert tag.weight == tag2.weight


def test_lha_weight():
    """ Test LHA weight tag. """
    tag = lhe.LHAweight()
    tag.identity = 0
    tag.add_weight('mur', 1)
    tag.add_weight('muf', 1)

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHAweight()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2
    assert tag.identity == tag2.identity
    assert tag.weights == tag2.weights


def test_lha_weight_group():
    """ Test LHA weight group tag. """
    tag = lhe.LHAweightgroup()
    tag.attrib['type'] = 'scale_variation'
    tag.attrib['combine'] = 'envelope'
    weights = [{'mur': 0.1E1, 'muf': 0.1E1},
               {'mur': 0.1E1, 'muf': 0.2E1},
               {'mur': 0.2E1, 'muf': 0.1E1},
               {'mur': 0.2E1, 'muf': 0.2E1},
               {'mur': 0.1E1, 'muf': 0.5E0},
               {'mur': 0.5E0, 'muf': 0.1E1},
               {'mur': 0.5E0, 'muf': 0.5E0}]
    tag.add_weights(weights)

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHAweightgroup()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2
    for key, _ in tag.children.items():
        assert tag.children[key] == tag2.children[key]


def test_lha_rwgt():
    """ Test LHA rwgt tag. """
    tag = lhe.LHArwgt()
    tag.add_weight(0, 0.783)

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHArwgt()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2
    for key, _ in tag.children.items():
        assert tag.children[key] == tag2.children[key]


def test_lha_init_rwgt():
    """ Test LHA initrwgt tag. """
    tag = lhe.LHAinitrwgt()

    group = lhe.LHAweightgroup()
    group.attrib['type'] = 'scale_variation'
    group.attrib['combine'] = 'envelope'
    weights = [{'mur': 0.1E1, 'muf': 0.1E1},
               {'mur': 0.1E1, 'muf': 0.2E1},
               {'mur': 0.2E1, 'muf': 0.1E1},
               {'mur': 0.2E1, 'muf': 0.2E1},
               {'mur': 0.1E1, 'muf': 0.5E0},
               {'mur': 0.5E0, 'muf': 0.1E1},
               {'mur': 0.5E0, 'muf': 0.5E0}]
    group.add_weights(weights)
    tag.add_weightgroup(group)

    weight = lhe.LHAweight()
    weight.identity = 0
    weight.add_weight('mur', 1)
    weight.add_weight('muf', 1)
    tag.add_weight(weight)

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHAinitrwgt()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2
    for key, _ in tag.children.items():
        assert tag.children[key] == tag2.children[key]


def test_lha_init():
    """ Test LHA init tag. """
    tag = lhe.LHAinit()
    tag.beams = [2212, 2212]
    tag.energy = [4000, 4000]
    tag.pdfgroup = [-1, -1]
    tag.pdfset = [21100, 21100]
    tag.wgtid = -4
    tag.text = 'This is a comment'
    process = lhe.LHAProcess(0.501090862e2, 0.89185414e-1, 0.50109093e2, 66)
    tag.processes.append(process)
    gen = lhe.LHAgenerator()
    gen.create('SomeGen1', '1.2.3', 'Some additional comments')
    tag.children['generator.SomeGen1'] = gen

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHAinit()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2
    for key, _ in tag.children.items():
        assert tag.children[key] == tag2.children[key]


def test_lha_event():
    """ Test LHA event tag. """
    tag = lhe.LHAevent()
    tag.attrib['npLO'] = ' -1 '
    tag.attrib['npNLO'] = ' 1 '
    tag.nprocess = 66
    tag.wgt = 0.501
    tag.scale = 0.14137
    tag.aqed = 0.755e-2
    tag.aqcd = 0.1211

    part = lhe.LHAParticle(5, -1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
    tag.particles.append(part)
    part = lhe.LHAParticle(2, -1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
    tag.particles.append(part)
    part = lhe.LHAParticle(24, 1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
    tag.particles.append(part)
    part = lhe.LHAParticle(2, 1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
    tag.particles.append(part)

    root = ElementTree.Element('root')
    tag.to_xml(root)

    xmlstr = ElementTree.tostring(root)

    root = ElementTree.fromstring(xmlstr)
    tag2 = lhe.LHAevent()
    head = root.findall(tag2.name)
    for child in head:
        assert child.tag == tag2.name
        tag2.from_xml(child)

    assert tag == tag2
    for key, _ in tag.children.items():
        assert tag.children[key] == tag2.children[key]


def test_lhef():
    """ Test LHEF file reader and writer. """
    import os

    test = lhe.LHEF('tmp.lhe')
    initrwgt = lhe.LHAinitrwgt()
    weightgroup = lhe.LHAweightgroup()
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

    init = lhe.LHAinit()
    init.beams = [2212, 2212]
    init.energy = [4000, 4000]
    init.pdfgroup = [-1, -1]
    init.pdfset = [21100, 21100]
    init.wgtid = -4
    init.text = 'this is a comment'
    process = lhe.LHAProcess(0.501090862e2, 0.89185414e-1, 0.50109093e2, 66)
    init.processes.append(process)
    gen = lhe.LHAgenerator()
    gen.create('SomeGen1', '1.2.3', 'Some additional comments')
    init.children['generator.SomeGen1'] = gen

    test.set_init_block(init)

    event = lhe.LHAevent()
    event.attrib['npLO'] = ' -1 '
    event.attrib['npNLO'] = ' 1 '
    event.nprocess = 66
    event.wgt = 0.501
    event.scale = 0.14137
    event.aqed = 0.755e-2
    event.aqcd = 0.1211

    part = lhe.LHAParticle(5, -1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
    event.particles.append(part)

    part = lhe.LHAParticle(2, -1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
    event.particles.append(part)

    part = lhe.LHAParticle(24, 1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
    event.particles.append(part)

    part = lhe.LHAParticle(2, 1, [0, 0, 0.1437E3, 1433E3, 0.48E1],
                           col=[501, 0], spin=0)
    event.particles.append(part)

    test.add_event(event)
    test.write()

    test2 = lhe.LHEF('tmp.lhe')
    test2.read()

    assert test.version == test2.version
    assert test.header.attrib == test2.header.attrib
    assert test.init == test2.init
    for i, _ in enumerate(test.events):
        assert test.events[i] == test2.events[i]
    assert test.header.text == test2.header.text

    os.remove('tmp.lhe')
