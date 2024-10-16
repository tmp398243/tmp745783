using Pkg: Pkg
using Lorenz63Filter
using Test
using TestReports
using TestReports.EzXML
using Aqua
using Documenter

function gen_runner_code(testsetname, testfilename, logfilename)
    runner_code = """
        using Pkg: Pkg

        Pkg.activate($(dirname(testfilename) |> repr))
        Pkg.develop(; path=$(joinpath(@__DIR__, "..") |> repr))
        Pkg.add(["Test", "TestReports"])
        Pkg.resolve()
        Pkg.instantiate()

        using Test: @testset
        using TestReports: TestReports, ReportingTestSet, any_problems, report
        using TestReports.EzXML: prettyprint

        ts = @testset ReportingTestSet $(repr(testsetname)) begin
            include($(repr(testfilename)))
        end

        # Flatten before calling `report` to avoid a `deepcopy`.
        flattened_testsets = TestReports.flatten_results!(ts)
        open($(repr(logfilename)), "w") do io
            prettyprint(io, report(flattened_testsets))
        end
        any_problems(flattened_testsets) && exit(TestReports.TESTS_FAILED)
        """
    return runner_code
end

errs = Vector{String}()
examples_dir = joinpath(@__DIR__, "..", "examples")

report_testsets = @testset ReportingTestSet "" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Lorenz63Filter; ambiguities=false)
        Aqua.test_ambiguities(Lorenz63Filter)
    end

    @info "Running package tests."
    @test true

    @info "Running doctests."
    # Set metadata for doctests.
    DocMeta.setdocmeta!(
        Lorenz63Filter, :DocTestSetup, :(using Lorenz63Filter, Test); recursive=true
    )
    doctest(Lorenz63Filter; manual=true)

    # Run examples.
    @info "Running examples."
    for example in readdir(examples_dir)
        @testset "Example: running $(example)" begin
            example_path = joinpath(examples_dir, example)
            script_path = joinpath(example_path, "main.jl")
            log_path = joinpath(example_path, "report.xml")
            tester_path = joinpath(mktempdir(), "tester.jl")
            @show tester_path example_path

            runner_code = gen_runner_code("Example: $(example)", script_path, log_path)
            open(tester_path, "w") do f
                return write(f, runner_code)
            end
            cmd = `$(Base.julia_cmd()) --color=no -- "$(tester_path)"`

            @info "Testing $example with \"$(tester_path)\""
            TestReports.runtests!(errs, example, cmd, log_path)
            @test length(errs) == 0 || errs[end] != example
        end
    end
end

xml_all = report(report_testsets)
a = xml_all
a_root = root(a)
a_attrs = Dict(c.name => parse(Int, c.content) for c in attributes(a_root))
for example in readdir(examples_dir)
    example_path = joinpath(examples_dir, example)
    log_path = joinpath(example_path, "report.xml")
    xml_ex = readxml(log_path)
    b = xml_ex

    b_root = root(b)
    b_elements = elements(b_root)
    a_nelems = length(elements(a_root))
    b_nelems = length(elements(b_root))
    for elem in eachelement(b_root)
        unlink!(elem)
        link!(a_root, elem)
    end
    @test length(elements(a_root)) == a_nelems + b_nelems

    for c in attributes(b_root)
        a_attrs[c.name] += parse(Int, c.content)
    end
end
for c in attributes(a_root)
    a_root[c.name] = a_attrs[c.name]
end

outputfilename = joinpath(@__DIR__, "..", "report.xml")
open(outputfilename, "w") do fh
    return print(fh, a_root)
end

exit(any_problems(report_testsets) || length(errs) > 0)
