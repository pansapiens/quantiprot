from distutils.core import setup
setup(
  name = 'quantiprot',
  version = '0.2.3',
  packages = ['quantiprot','quantiprot.utils','quantiprot.metrics','quantiprot.analysis'],
  package_data = { 'quantiprot.metrics': [ 'data/aaindex1'] },
  description = 'Quantitative characteristics of protein sequences',
  author = 'Witold Dyrka, Bogumil M. Konopka, Marta Marciniak',
  author_email = 'witold.dyrka@pwr.edu.pl',
  url = 'https://git.e-science.pl/wdyrka/quantiprot',
  download_url = 'https://git.e-science.pl/wdyrka/quantiprot/repository/archive.tar.gz?ref=0.2.3',
  keywords = ['protein', 'sequence', 'quantitative', 'analysis', 'rqa', 'n-gram'],
  classifiers = [],
  install_requires = ['numpy>=1.11.0', 'requests>=2.10.0'],
)
