
import PALEOboxes as PB

model =  PB.create_model_from_config("configfluxes.yaml", "model1")

PB.TestUtils.bench_model(model)