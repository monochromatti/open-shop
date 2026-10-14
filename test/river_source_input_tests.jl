@testset "Direct source compilation and input domains" begin
    f=mixed_source_small_fixture();c=f.case;s=c.system
    @test case_dict(OpenSHOP._river_replace(c;system=normalize_river_connections(s)))==case_dict(c)
    badplant=OpenSHOP._river_replace(s.plants[1];discharge_river=:Missing)
    badsystem=OpenSHOP._river_replace(s;plants=vcat([badplant],s.plants[2:end]))
    @test_throws ArgumentError normalize_river_connections(badsystem)
    @test_throws ArgumentError OpenSHOP._river_replace(c;system=badsystem)
    # Builders also validate manually assembled Julia systems. A redirected
    # discharge cannot disappear into a reach whose source is a reservoir.
    badplant=OpenSHOP._river_replace(s.plants[1];discharge_river=:ReceiverOutlet)
    badsystem=OpenSHOP._river_replace(s;plants=vcat([badplant],s.plants[2:end]))
    @test_throws ArgumentError OpenSHOP.validate_inputs(OpenSHOP._river_replace(c;system=badsystem))
    public=case_dict(c)
    wrong=deepcopy(public);wrong["rivers"][1]["source"]=:UpperA
    @test_throws ArgumentError case_from_dict(wrong)
    wrong=deepcopy(public);wrong["rivers"][1]["law"]=:controlled
    @test_throws ArgumentError case_from_dict(wrong)
    wrong=deepcopy(public);wrong["rivers"][1]["inflow"]=-0.1
    @test_throws ErrorException case_from_dict(wrong)
    wrong=deepcopy(public);wrong["rivers"][1]["inflow"]=9.0
    @test_throws ArgumentError case_from_dict(wrong)

    # A natural source need not have a static positive default: its operating
    # series supplies water, and the zero-storage source is still internal.
    river=River(name=:Creek,target=:Sea,capacity=5.0,water_value=0.0,
        deterministic_delay=0.25,curves=RiverRouting.DelayCurve[],
        history_grid=[-1.0,0.0],history_release=[0.0])
    system=normalize_river_connections(HydroSystem(reservoirs=Reservoir[],
        junctions=Junction[],boundaries=[Boundary(name=:Sea,head=0.0)],
        tunnels=Tunnel[],plants=Plant[],generators=Generator[],rivers=[river]))
    natural=ScheduleCase(name="operational natural inflow",system=system,
        grid=[0.0,1.0,2.0],prices=zeros(2),operations=[OperationalSeries(
            object=:Creek,attribute=:inflow,times=[0.0,1.0],values=[2.0,4.0])])
    @test OpenSHOP.validate_inputs(natural)
    z=simulate(natural,zeros(0,2),zeros(1,2))
    @test z["converged"]
    @test z["river_release"]≈[2.0 4.0] atol=1e-12
    @test z["arrival_volume"]≈0.0036 .* [1.5 3.5] atol=1e-12
    @test z["transit"]≈0.0036 .* [0.0 0.5 1.0] atol=1e-12
    @test isempty(case_dict(natural)["river_junctions"])
    cached=deterministic_network_data(natural)
    stopped=OpenSHOP._river_replace(natural;operations=[OperationalSeries(
        object=:Creek,attribute=:inflow,times=[0.0],values=[0.0])])
    @test_throws ArgumentError OpenSHOP._transport_data(stopped,cached)

end
