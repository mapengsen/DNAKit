;;; GNU Guix package definition for DNAKit.

(use-modules (gnu packages check)
             (gnu packages python-build)
             (gnu packages python-xyz)
             (guix build-system pyproject)
             (guix gexp)
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
  (arguments
   (list
    ;; Tests import the source-only benchmarks and validation modules.
    #:test-flags #~(list "-o" "pythonpath=.")
    #:phases
    #~(modify-phases %standard-phases
        (add-after 'unpack 'patch-test-shells
          (lambda _
            ;; The test executables are generated after patch-shebangs.
            (substitute* '("tests/unit/test_external_backend_adapters.py"
                           "tests/unit/test_remote_download.py")
              (("#!/bin/sh")
               (string-append "#!" (which "sh")))))))))
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
