using RegressionTests, Chairmarks

load_time = @elapsed using OffsetTables
@track load_time

x = rand(10)
ttf_construction = @elapsed ot = OffsetTable(x)
@track ttf_construction

ttf_sample = @elapsed rand(ot)
@track ttf_sample

ttf_sample_multi = @elapsed rand(ot, 30)
@track ttf_sample_multi

for n in [1, 10, 100, 1000]
    @track @b rand(n) OffsetTable seconds=.02
    @track @b OffsetTable(rand(n)) rand seconds=.02
    @track @b OffsetTable(rand(n)) rand(_, 100) seconds=.02
end

for n in [3, 30, 300]
    @track @b OffsetTable(rand(n)) hash seconds=.01
    @track @b OffsetTable(rand(n)), OffsetTable(rand(n)) (==)(_...) seconds=.01
    @track @b (x = OffsetTable(rand(n)); (x, deepcopy(x))) (==)(_...) seconds=.01
end
