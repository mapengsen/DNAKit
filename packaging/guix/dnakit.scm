;;; GNU Guix package definition for DNAKit.

(use-modules (gnu packages check)
             (gnu packages python-build)
             (gnu packages python-xyz)
             (guix build-system pyproject)
             (guix download)
             ((guix licenses) #:prefix license:)
             (guix packages))

(package
  (name "python-dnakit")
  (version "0.1.1")
  (source
   (origin
     (method url-fetch)
     (uri (pypi-uri "dnakit" version))
     (sha256
      (base32 "1jyz4fiibh21nwq2gdi4i4xjz5ryii9bnvd3xgxrzmj9kc0idmp7"))))
  (build-system pyproject-build-system)
  (native-inputs
   (list python-numpy
         python-pytest
         python-setuptools
         python-wheel))
  (propagated-inputs
   (list python-pyyaml
         python-rich
         python-tomli
         python-typer))
  (home-page "https://github.com/mapengsen/DNAKit")
  (synopsis "Deterministic tools for DNA sequence analysis")
  (description
   "DNAKit is a reproducible Python toolkit for DNA sequence normalization,
file conversion, descriptors, pattern scanning, similarity analysis, dataset
preparation, evaluation, simulation, and visualization.")
  (license license:expat))
