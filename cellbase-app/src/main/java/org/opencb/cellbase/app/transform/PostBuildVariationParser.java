package org.opencb.cellbase.app.transform;

import com.beust.jcommander.ParameterException;
import com.fasterxml.jackson.databind.MapperFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectReader;
import com.fasterxml.jackson.databind.ObjectWriter;
import org.apache.commons.collections.map.HashedMap;
import org.opencb.biodata.formats.variant.io.JsonVariantReader;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.biodata.models.variant.avro.PopulationFrequency;
import org.opencb.biodata.models.variant.avro.VariantAnnotation;
import org.opencb.cellbase.core.serializer.CellBaseSerializer;
import org.rocksdb.Options;
import org.rocksdb.RocksDB;
import org.rocksdb.RocksDBException;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

/**
 * Created by fjlopez on 08/06/16.
 */
public class PostBuildVariationParser extends CellBaseParser {

    private final Path variationPath;
    private final Path annotationPath;
    private final Path frequenciesPath;

    private static final int BATCH_SIZE = 1000;
    private static final String ANNOTATION_FIELD = "annotation";
    private static final String FREQUENCIES_FIELD = "annotation.populationFrequencies";

    private final Map<String, ObjectReader> byteReaderMap = new HashedMap(2);

    public PostBuildVariationParser(Path variationPath, Path annotationPath, Path frequenciesPath,
                                     CellBaseSerializer serializer) {
        super(serializer);
        this.variationPath = variationPath;
        this.annotationPath = annotationPath;
        this.frequenciesPath = frequenciesPath;
        ObjectMapper mapper = new ObjectMapper();
        byteReaderMap.put(ANNOTATION_FIELD,
                mapper.readerFor(mapper.getTypeFactory().constructType(VariantAnnotation.class)));
        byteReaderMap.put(FREQUENCIES_FIELD,
                mapper.readerFor(mapper.getTypeFactory().constructCollectionType(List.class, PopulationFrequency.class)));
    }

    public void parse() throws IOException {
        RocksDB annotationIdx = indexVariantField(annotationPath, ANNOTATION_FIELD);
        RocksDB frequenciesIdx = indexVariantField(frequenciesPath, FREQUENCIES_FIELD);

        JsonVariantReader variantReader = new JsonVariantReader(variationPath.toString());
        variantReader.open();
        List<Variant> variantList = variantReader.read(BATCH_SIZE);
        while (!variantList.isEmpty()) {
            for (Variant variant : variantList) {
                if (annotationIdx != null) {
                    mergeAnnotation(variant.getAnnotation(), getIndexedValue(annotationIdx,
                            byteReaderMap.get(ANNOTATION_FIELD), variant.toString()));
                }
                if (frequenciesIdx != null) {
                    variant.getAnnotation().setPopulationFrequencies(getIndexedValue(frequenciesIdx,
                            byteReaderMap.get(FREQUENCIES_FIELD), variant.toString()));
                }
                serializer.serialize(variant);
            }
            variantList = variantReader.read(BATCH_SIZE);
        }
        annotationIdx.close();
        frequenciesIdx.close();
        variantReader.close();
        serializer.close();
    }

    private void mergeAnnotation(VariantAnnotation destination, VariantAnnotation origin) {
        if (origin != null) {
            destination.setChromosome(origin.getChromosome());
            destination.setStart(origin.getStart());
            destination.setReference(origin.getReference());
            destination.setAlternate(origin.getAlternate());
            destination.setDisplayConsequenceType(origin.getDisplayConsequenceType());
            destination.setConsequenceTypes(origin.getConsequenceTypes());
            destination.setConservation(origin.getConservation());
            destination.setGeneExpression(origin.getGeneExpression());
            destination.setGeneTraitAssociation(origin.getGeneTraitAssociation());
            destination.setGeneDrugInteraction(origin.getGeneDrugInteraction());
            destination.setVariantTraitAssociation(origin.getVariantTraitAssociation());
            destination.setFunctionalScore(origin.getFunctionalScore());
        }
    }

//    private <T> T getIndexedValue (RocksDB db, String key, Class<T> clazz) {
//        return getIndexedValue(db, key, );
//    }
//    private <T> List<T> getIndexedValueList (RocksDB db, String key, Class<T> clazz) {
////        JavaType type = mapper.getTypeFactory()
////                .constructParametrizedType(List.class, List.class, clazz);
//
//        JavaType type = mapper.getTypeFactory()
//                .constructCollectionType(List.class, clazz);
//        return getIndexedValue(db, key, type);
//    }

    private <T> T getIndexedValue(RocksDB db, ObjectReader byteReader, String key) {
        byte[] dbContent;
        try {
            dbContent = db.get(key.getBytes());
            if (dbContent == null) {
                return null;
            } else {

                return byteReader.readValue(dbContent);
            }
        } catch (RocksDBException | IOException e) {
            e.printStackTrace();
            return null;
        }
    }

    private RocksDB indexVariantField(Path path, String field) {
        Object[] dbConnection = getDBConnection(path.toString() + ".idx");
        RocksDB db = (RocksDB) dbConnection[0];
        boolean indexingNeeded = (boolean) dbConnection[1];
        if (indexingNeeded) {
            logger.info("Creating index DB at {}.idx", path.toString());
            ObjectMapper jsonObjectMapper = new ObjectMapper();
            jsonObjectMapper.configure(MapperFeature.REQUIRE_SETTERS_FOR_GETTERS, true);
            ObjectWriter jsonObjectWriter = jsonObjectMapper.writer();
            try {
                JsonVariantReader variantReader = new JsonVariantReader(path.toString());
                variantReader.open();
                List<Variant> variantList = variantReader.read(BATCH_SIZE);
                int variantCounter = 0;
                while (!variantList.isEmpty()) {
                    for (Variant variant : variantList) {
                        byte[] toIndex;
                        switch (field) {
                            case ANNOTATION_FIELD:
                                toIndex = jsonObjectWriter.writeValueAsBytes(variant.getAnnotation());
//                                toIndex = variant.getAnnotation().toString().getBytes();
                                break;
                            case FREQUENCIES_FIELD:
                                toIndex = jsonObjectWriter.writeValueAsBytes(variant.getAnnotation()
                                        .getPopulationFrequencies());
                                break;
                            default:
                                throw new ParameterException("No action defined for indexing field: " + field);
                        }
                        db.put((variant.toString()).getBytes(), toIndex);
                        variantCounter++;
                        if (variantCounter % 100000 == 0) {
                            logger.info("{} variants indexed", variantCounter);
                        }
                    }
                    variantList = variantReader.read(BATCH_SIZE);
                }
                variantReader.close();
            } catch (IOException | RocksDBException e) {
                e.printStackTrace();
                System.exit(1);
            }
        } else {
            logger.info("Index found at {}.idx", annotationPath);
            logger.info("Skipping index creation");
        }

        return db;
    }

    private Object[] getDBConnection(String dbLocation) {
        boolean indexingNeeded = !Files.exists(Paths.get(dbLocation));
        // a static method that loads the RocksDB C++ library.
        RocksDB.loadLibrary();
        // the Options class contains a set of configurable DB options
        // that determines the behavior of a database.
        Options options = new Options().setCreateIfMissing(true);
        RocksDB db = null;
        try {
            // a factory method that returns a RocksDB instance
            if (indexingNeeded) {
                db = RocksDB.open(options, dbLocation);
            } else {
                db = RocksDB.openReadOnly(options, dbLocation);
            }
            // do something
        } catch (RocksDBException e) {
            // do some error handling
            e.printStackTrace();
            System.exit(1);
        }

        return new Object[]{db, indexingNeeded};

    }

}
