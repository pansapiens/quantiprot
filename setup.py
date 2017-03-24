from distutils.core import setup
setup(
  name = 'quantiprot',
  packages = ['quantiprot','quantiprot.utils','quantiprot.metrics','quantiprot.analysis'],
  version = '0.2.2',
  description = 'Quantitative characteristics of protein sequences',
  author = 'Witold Dyrka, Bogumil M. Konopka, Marta Marciniak',
  author_email = 'witold.dyrka@pwr.edu.pl',
  url = 'https://git.e-science.pl/wdyrka/quantiprot',
  download_url = 'https://git.e-science.pl/wdyrka/quantiprot/repository/archive.tar.gz?ref=0.2.2',
  keywords = ['protein', 'sequence', 'quantitative', 'analysis', 'rqa', 'n-gram'],
  classifiers = [],
  install_requires = ['numpy>=1.11.0', 'requests>=2.10.0'],
)
