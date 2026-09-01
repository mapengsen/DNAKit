;;; GNU Guix package definition for DNAKit.

(use-modules (gnu packages check)
             (gnu packages python-build)
             (gnu packages python-xyz)
             (guix build-system pyproject)
             (guix git-download)
             ((guix licenses) #:prefix license:)
             (guix packages))

(package
  (name "python-dnakit")
  (version "0.1.1")
  (source
   (origin
     (method git-fetch)
     (uri (git-reference
           (url "https://github.com/mapengsen/DNAKit")
           (commit (string-append "v" version))))
     (file-name (git-file-name name version))
     (sha256
      (base32 "1cra4r62x0fzwsjlbhc3gl5985dlh3p74y9v77idbmbnd76jmv2w"))))
  (build-system pyproject-build-system)
  (native-inputs
   (list python-numpy
         python-pytest
         python-setuptools
         python-wheel))
  (propagated-inputs
   (list python-pyyaml
         python-rich
         python-typer))
  (home-page "https://github.com/mapengsen/DNAKit")
  (synopsis "Deterministic tools for DNA sequence analysis")
  (description
   "DNAKit is a reproducible Python toolkit for DNA sequence normalization,
file conversion, descriptors, pattern scanning, similarity analysis, dataset
preparation, evaluation, simulation, and visualization.")
  (license license:expat))
