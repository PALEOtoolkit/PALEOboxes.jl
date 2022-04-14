import PALEOboxes as PB

include("ReactionRateStoichMock.jl")


model = PB.create_model_from_config("configratestoich.yaml", "model1")

PB.TestUtils.bench_model(model)
