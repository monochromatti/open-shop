@testset "A numerical mesh projection does not execute invalid controls" begin
    c=warm_dispatch_fixture()
    u=[1 1;0 0]
    x=dispatch_from_controls(c,u,[4.0 4.0;0.0 0.0],ones(1,2))
    guess=deepcopy(x)
    guess.generator_q[1,:].=1e6
    grid=[0.0,0.25,1.0,1.6,2.0]
    projected,_=OpenSHOP._dispatch_start(c,u,guess;grid)
    @test projected!==nothing
    @test projected.storage[:,[1,end]]≈guess.storage[:,[1,end]]
    @test sum(projected.arrival_volume;dims=2)≈sum(guess.arrival_volume;dims=2)
    @test projected.power*diff(grid)≈guess.power*diff(guess.grid)
    @test all(guess.generator_q[1,:].==1e6)
    restored=solve_verified(c;u,warm=guess,feasibility_only=true,
        time_limit=10.0,deadline=time()+30.0)
    @test restored.accepted
    if restored.accepted
        @test validate(c,restored.solution)["valid"]
        @test replay_audit(c,restored.solution)["valid"]
        @test maximum(restored.solution.generator_q)<=5.0+1e-6
    end
end

@testset "An interrupted integration cannot be accepted as complete" begin
    c=warm_dispatch_fixture()
    u=[1 1;0 0]
    q=[4.0 4.0;0.0 0.0]
    gate=ones(1,2)
    stopped=simulate(c,q,gate;commitment=u,deadline=time()-1)
    @test !stopped.converged
    @test_throws ArgumentError dispatch_from_controls(c,u,q,gate;deadline=time()-1)
    @test simulate(c,q,gate;commitment=u).converged
    complete=dispatch_from_controls(c,u,q,gate)
    @test !replay_audit(c,complete;deadline=time()-1)["valid"]
end
