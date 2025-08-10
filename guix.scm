;; guix.scm for povu
;;
;; To use this file to build HEAD of povu:
;;
;;   guix build -f guix.scm
;;
;; To get a development container:
;;
;;   guix shell -C -D -f guix.scm
;;
;; and build:
;;
;;   rm -rf build
;;   mkdir build && cd build
;;   cmake -DCMAKE_BUILD_TYPE=Debug ..
;;   cmake --build . --verbose -- -j $(nproc)
;;
;; For release build:
;;
;;   cmake -DCMAKE_BUILD_TYPE=Release ..
;;   cmake --build . --verbose -- -j $(nproc)
;;
;; For static build:
;;
;;   cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC=ON ..
;;   cmake --build . --verbose -- -j $(nproc)

;; by Uncle Claude & Andrea Guarracino (c) 2025

(use-modules
  (ice-9 popen)
  (ice-9 rdelim)
  ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix download)
  (guix git-download)
  (guix build-system cmake)
  (guix utils)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages build-tools)
  (gnu packages cmake)
  (gnu packages commencement)
  (gnu packages cpp)
  (gnu packages gcc)
  (gnu packages gdb)
  (gnu packages pkg-config)
  (gnu packages pretty-print)
  (gnu packages version-control)
)

(define %source-dir (dirname (current-filename)))

(define %git-commit
  (read-string (open-pipe "git show HEAD 2>/dev/null | head -1 | cut -d ' ' -f 2 || echo unknown" OPEN_READ)))

(define-public povu-git
  (package
    (name "povu-git")
    (version (git-version "0.1.0" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (inputs
     `(("coreutils" ,coreutils)
       ("fmt" ,fmt)  ; Alternative to std::format for better compatibility
       ("gcc" ,gcc-13)
       ("gcc-lib" ,gcc-13 "lib")
       ("gcc-toolchain" ,gcc-toolchain)
       ("gdb" ,gdb)
       ("git" ,git)
       ("zlib" ,zlib)
       ("bzip2" ,bzip2)
       ("xz" ,xz)
       ))
    (native-inputs
     `(("cmake" ,cmake)
       ("pkg-config" ,pkg-config)
       ))
    (arguments
      `(#:configure-flags
        (list "-DCMAKE_BUILD_TYPE=Release")
        #:tests? #f  ; Disable tests for now
        #:phases
        (modify-phases
         %standard-phases
         ;; Store build version in the executable
         (add-after 'unpack 'set-version
           (lambda _
             (format #t "Building povu version ~a~%" ,version)
             #t))
         )
        ))
     (synopsis "povu - pangenome graph variation toolkit")
     (description
"povu is a lightweight tool for exploring regions of genomic variation
in pangenome graphs. It provides algorithms for graph manipulation,
alignment, and variant calling from pangenome representations.")
     (home-page "https://github.com/pangenome/povu")
     (license license:expat)))

povu-git