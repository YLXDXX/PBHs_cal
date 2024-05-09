add_rules("mode.debug", "mode.release")

target("./PBHs") --Name of the generated executable file
    set_kind("binary")
    add_files("src/main.c")
    add_files("src/func/*.c")
    add_files("src/func/Math_Method/*.c")
    add_files("src/func/General/*.c")
    add_files("src/func/Peak_Theory/*.c")
    add_files("src/func/PS_Method/*.c")
    add_files("src/func/Other/*.c")
    add_files("src/func/GW_induced/*.c")
    add_files("src/func/Test/*.c")
    add_files("src/set/*.c")

--
-- If you want to known more usage about xmake, please see https://xmake.io
--
-- ## FAQ
--
-- You can enter the project directory firstly before building project.
--
--   $ cd projectdir
--
-- 1. How to build project?
--
--   $ xmake
--
-- 2. How to configure project?
--
--   $ xmake f -p [macosx|linux|iphoneos ..] -a [x86_64|i386|arm64 ..] -m [debug|release]
--
-- 3. Where is the build output directory?
--
--   The default output directory is `./build` and you can configure the output directory.
--
--   $ xmake f -o outputdir
--   $ xmake
--
-- 4. How to run and debug target after building project?
--
--   $ xmake run [targetname]
--   $ xmake run -d [targetname]
--
-- 5. How to install target to the system directory or other output directory?
--
--   $ xmake install
--   $ xmake install -o installdir
--
-- 6. Add some frequently-used compilation flags in xmake.lua
--
-- @code
--    -- add debug and release modes
--    add_rules("mode.debug", "mode.release")
--
--    -- add macro definition
--    add_defines("NDEBUG", "_GNU_SOURCE=1")
--
--    -- set warning all as error
    set_warnings("all", "error")
--
--    -- set language: c99, c++11
--    set_languages("c99", "c++11")
--
--    -- set optimization: none, faster, fastest, smallest
    set_optimize("fastest")
--
--    -- add include search directories
    add_includedirs("/usr/local/include/flint")
--
--    -- add link libraries and search directories
--    add_links("tbox")
--    add_linkdirs("/usr/local/lib", "/usr/lib")
--
--    -- add system link libraries
--    add_syslinks("z", "pthread")
--
--    -- add compilation and link flags
    add_cxflags("-Wall")
    add_ldflags("-Wall", "-lflint", "-lm",{force = true})
--
-- @endcode
--

