((nil . ((projectile-project-compilation-cmd . "cmake -Bbuild -DPOVU_ENABLE_TESTING=ON -DCMAKE_BUILD_TYPE=Debug -S. && cmake --build build")
         (projectile-project-run-cmd . "./bin/povu")
         (projectile-project-test-cmd . "ctest --verbose --test-dir build --output-on-failure"))))
