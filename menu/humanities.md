---
layout: page
title: Humanities Dataset Details
---

**You will be able to download the dataset on October 1, 2018.**

The Orlando British Women's Writing Dataset provides a rich set of linked open data representing women's
literary history from the beginnings to the present, concentrated on writing in English in the British
Isles but with tentacles out to other languages, literary traditions, and parts of the world. It emerges
from the ongoing experiments in literary history being conducted by
[the Orlando Project](http://www.artsrn.ualberta.ca/orlando), whose textbase is published and regularly
updated and augmented as
[Orlando: Women’s Writing in the British Isles from the Beginnings to the Present](http://orlando.cambridge.org)
by Cambridge University Press since 2006, and from the
[Canadian Writing Research Collaboratory](https://beta.cwrc.ca)'s work in Linked Open Data.

The Orlando Textbase is a semi-structured collection of biocritical entries providing detailed
information on the lives and writing of more than 1400 writers with accompanying literary, social, and
political materials to provide context to its representation of literary history. It does not contain
digitized versions of primary texts.

The internal linking of biographical information is using the
[Canadian Writing Research Collaboratory (CWRC) ontology](http://sparql.cwrc.ca/ontology/cwrc.html), with
selective linking out to other ontology terms and linked data entities.

## How to read the data

The dataset contains more than 2 million triples stored in the
[Resource Description Framework](https://www.w3.org/RDF) (RDF) format. RDF is a standard format for
storing hierarchical or linked data called *graphs* in RDF.

### Learning RDF triplestore

For a good introduction to RDF for those without any prior experience, we recommend saving an example
from
[this Mozilla tutorial](https://developer.mozilla.org/en-US/docs/Mozilla/Tech/XUL/Tutorial/Introduction_to_RDF)
as an RDF file `zoo.rdf` and then reading and analyzing it from Python using the
[RDFLib library](https://github.com/RDFLib/rdflib). This RDF results in 13 triples describing the
relationship between objects and their properties.

```python
from rdflib import Graph, RDF, URIRef, RDFS
g = Graph()   # create an empty graph to load data into
result = g.parse("zoo.rdf")   # load the data
print("graph has %s statements" % len(g))   # 13 statements
```

You can then print all the triples:

```python
for subject, predicate, object in g:
    print(subject, predicate, object)
```

resulting in the following output (after manual cleanup and reordering):

```bash
all-animals is a class for ordered components
all-animals for element 1 has lion
all-animals for element 2 has tarantula
all-animals for element 3 has hippopotamus
lion is named Lion
lion is of species Panthera leo
lion is of class Mammal
tarantula is named Tarantula
tarantula is of species Avicularia avicularia
tarantula is of class Arachnid
hippopotamus is named Hippopotamus
hippopotamus is of species Hippopotamus amphibius
hippopotamus is of class Mammal
```

These triples completely describes all information stored in the RDF-XML file, in a human-readable
form. You can also list all predicates of a given subject

```python
predicates = g.predicates(subject=URIRef('http://www.some-fictitious-zoo.com/mammals/lion'))
for j, p in enumerate(predicates):
    print(j, p)
```

or print all unique subjects, predicates, objects

```python
subjects = set(g.subjects())
for s in subjects:
    print(s)

predicates = set(g.predicates())
for p in predicates:
    print(p)

objects = set(g.objects())
for o in objects:
    print(o)
```

and then further filter these data with Python. The 13 triples completely describe all data objects and
their connections in the file, so that one could use this information for plotting or reconstructing the
original RDF. Note that the source RDF file takes much less space then the resulting triples.

### Competition data

The competition dataset contains 5475 files across three directories

```bash
Bibliography/ Biography/ CulturalForms/
```

Some of the files are stored in both RDF-XML and RDF-TTL formats, so you can choose the format you
like. In addition to using separate generators `g.subjects()`, `g.predicates()`, `g.objects()` as
outlined above, you can also list triples that match a given pattern, e.g., a given subject with all
possible predicates and objects:

```python
# list all predicates and objects of a given subject
for s, p, o in g.triples((URIRef("http://cwrc.ca/cwrcdata/Abdy_Maria"), None, None)):
    print(p, o)
```

The bibliography file `Bibliography.ttl` is actually a RDF-TTL file, and you can read it with
RDFLib specifying the format explicitly

```python
from rdflib import Graph
g = Graph()   # create an empty graph to load data into
result = g.parse("Bibliography/Bibliography.ttl", format="ttl")
```

For the list of variables, their description, and how they interact with each other, it is best to
consult the [CWRC Ontology Specification](http://sparql.cwrc.ca/ontology/cwrc.html) and the
[CWRC preamble](http://sparql.cwrc.ca/ontologies/cwrc-preamble-EN.html). The dataset also adops some
external namespaces, classes and terms listed
[here](http://sparql.cwrc.ca/ontologies/cwrc-preamble-EN.html#linkages). At the beginning of each RDF
file you will find all the external vocabularies that are used in that file, e.g.

```bash
$ head Biography/RDF-XML/abdyma.rdf
<?xml version="1.0" encoding="utf-8"?>
<rdf:RDF
  xmlns:dcterms="http://purl.org/dc/terms/"
  xmlns:foaf="http://xmlns.com/foaf/0.1/"
  xmlns:oa="http://www.w3.org/ns/oa#"
  xmlns:org="http://www.w3.org/ns/org#"
  xmlns:cwrc="http://sparql.cwrc.ca/ontologies/cwrc#"
  xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
>
```

Consulting these sources in relation to the data is key to understanding the dataset.

## Research questions

To demonstrate the broad type of information in the dataset, the CWRC ontology provides a number of
[competency questions](http://sparql.cwrc.ca/ontologies/cwrc-preamble-EN.html#competency). For the
purpose of this competition, we suggest some more focused questions that the current dataset is
positioned to illuminate:

- What is the relevance of family size to religious affiliation, socioeconomic status, or other factors?
- Does social identity become more diverse over time? Are hybrid cultural forms more common than they
  were previously? Where are the clusters of hybridity and what are the outliers?
- What do bibliometric visualizations reveal about the publishing networks in the bibliographical data?
- Does the rise and fall of various genres over time within the works of Orlando authors correlate with
  scholarly accounts of when those genres rose and fell?
- Which authors are the most isolated in the network, and what factors (for instance, social class, place
  of publication) seem to be associated with being an outlier?
- Can number of publications be related to the number of children a woman writer had?
- Who are the most connected authors in the textbase (most points of contact with other authors) and what
  are the most common types of connections, in terms of either specific relationships or the contexts in
  which they occur? Do the types of connections tracked in the dataset modify over time?
- What can a visualization reveal both about the structure of the ontology and the data associated with
  different components of it? For instance, can it reflect the number and types of instances associated
  with different components of the ontology, such as by showing the shifting proportions of different
  Biographical Contexts such as religion, politics, and sexuality as they occur in entries over time?

These questions could be approached by visualizations in a range of formats including charts, network
graphs, geospatial maps, heat maps, trend visualization, and infographics.

#### An example of a chart approach

For the research question "What is the relevance of family size to religious affiliation, socioeconomic
status, or other factors?", produce a series of charts (bar, line, pie charts, etc.) that answer this
question. These could be 2D or 3D. Make the process interactive so that the user can select the factors
for a multi-dimensional visualization.

#### An example of a network graph approach

What kinds of biographical networks connect British women writers to each other and to other writers? How
extensive are kinship networks as opposed to networks based on political or religious affiliation?

An example of a map approach: Can you map out all the geographic information relevant to a person in
Orlando? Can you also map out the person and all the people that they are connected to in Orlando? Does
this map information change over time period, with ethnicity, religion, etc.? How do the countries
associated with geographical heritage change over time?

#### An example of a heat map approach

Taking one of the research questions provided above such as "Does social identity become more diverse
over time?", can you represent this as a series of heat maps that plot out social identity (religion,
ethnicity, etc.) over time?  This could even be in the form of a video such as
[this one](https://www.flickr.com/photos/150411108@N06/43350961005/#).

#### An example of trend analysis

The research question "How is the number of publications related to the number of children a woman writer
had?" could be illustrated using dynamic graphs, examples of which can be seen in videos of Hans
Rosling's Gapminder visualizations such as [this one](https://www.youtube.com/watch?v=jbkSRLYSojo).

#### An example of an infographics approach

If you Google Emily Brontë, you will likely come across this infographic:

![alt text]({{ site.baseurl }}/assets/img/bronte.png "Emily Brontë")

Using the Orlando data, can you produce a more detailed infographic of the author? Can you then expand
this to have additional infographics of people who are connected to her and who also appear in the
Orlando data? Can you find information in the Orlando data and the ontologies that would produce a much
different infographic than Google's? What would an infographic of a larger group, such as all writers
from a particular historical period, or all writers of a given genre of literature, look like?











## Sample Visualizations

We will provide sample visualizations shortly.






<!-- To give you an idea of the type of data in these files and to help you with actual visualizations, we -->
<!-- provide two sample visualizations, one done with ParaView and the other one with VisIt. Both workflows -->
<!-- demonstrate loading of all 11 VTK files. -->

<!-- The ParaView state file bladesWithLines.pvsm stores the pipeline to visualize the blades (coloured by the -->
<!-- pressure on their surfaces) and the airflow around them with uniform-colour streamlines. You can point -->
<!-- ParaView to this state file with File - Load State..., or start ParaView from the command line with -->
<!-- "paraview --state=bladesWithLines.pvsm". The resulting image bladesWithLines.png is shown below. -->

<!-- /files/webfm/Communications/bladesWithLines.png -->

<!-- The VisIt Python script positiveNegativePressure.py renders semi-transparent isosurfaces of positive -->
<!-- (blue) and negative (turquoise) pressure around the blades. You can run this script in VisIt either from -->
<!-- Controls - Launch CLI... or from Controls - Command..., or from the command line with "visit -nowin -cli -->
<!-- -s positiveNegativePressure.py". The resulting image positiveNegativePressure0000.png is shown below. -->

<!-- /files/webfm/Communications/positiveNegativePressure0000.png -->

<!-- We are looking for innovative visualizations of this dataset. For example, one could enhance these -->
<!-- renderings by drawing streamlines around the isosurfaces and producing some animations such as spinning -->
<!-- the visualization around the vertical axis or gradually turning on/off various visualization -->
<!-- elements. Speaking more generally, a nice animation would help us explore the spatial range and values of -->
<!-- multiple variables and show how various elements of the simulation are tied together. -->
