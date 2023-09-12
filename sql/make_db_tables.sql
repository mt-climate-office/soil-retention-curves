-- soil.model_error definition

-- Drop table

-- DROP TABLE soil.model_error;

CREATE TABLE soil.model_error (
	station text NOT NULL,
	"depth" numeric NOT NULL,
	model text NOT NULL,
	mse numeric NOT NULL
);
CREATE UNIQUE INDEX soil_error_idx ON soil.model_error USING btree (station, depth, model);

-- soil.parameters definition

-- Drop table

-- DROP TABLE soil.parameters;

CREATE TABLE soil.parameters (
	station text NOT NULL,
	"depth" numeric NOT NULL,
	model text NOT NULL,
	params jsonb NOT NULL,
	CONSTRAINT soil_params_unique PRIMARY KEY (station, depth, model)
);
CREATE UNIQUE INDEX soil_params_idx ON soil.parameters USING btree (station, depth, model);

-- soil.raw_soil_data definition

-- Drop table

-- DROP TABLE soil.raw_soil_data;

CREATE TABLE soil.raw_soil_data (
	station text NOT NULL,
	"depth" numeric NOT NULL,
	kpa numeric NOT NULL,
	vwc numeric NOT NULL
);
CREATE INDEX soil_raw_idx ON soil.raw_soil_data USING btree (station, dept