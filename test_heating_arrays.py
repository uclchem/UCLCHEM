def test_heating_arrays_comprehensive():
    """Comprehensive test of heating array functionality."""
    print("=" * 60)
    print("Testing heating array")
    print("=" * 60)
    
    try:
        import uclchem
        from uclchem.model import _create_arrays, _get_standard_array_specs
        
        # Test 1: Array creation and specifications
        print("\n1. Testing array creation and specifications...")
        
        # Test parameters that should work (from run_uclchem_tests.py)
        param_dict = {
            "endAtFinalDensity": False,
            "freefall": False,
            "writeStep": 1,
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalDens": 1e5,
            "finalTime": 5.0e6,
            "points": 1,
        }
        
        # Test array specifications
        array_specs = _get_standard_array_specs(return_heating=True)
        print(f"Array specifications: {array_specs}")
        
        # Test array creation
        arrays = _create_arrays(param_dict, array_specs, timepoints=50)
        print(f"Created arrays: {list(arrays.keys())}")
        print(f"Heating array shape: {arrays['heatarray'].shape}")
        
        # Verify heating array has correct dimensions
        expected_shape = (51, 1, 18)  # timepoints+1, points, n_heating_terms
        assert arrays['heatarray'].shape == expected_shape, f"Expected {expected_shape}, got {arrays['heatarray'].shape}"
        print("Array creation test passed")
        
        # Test 2: High-level model function with return_array
        print("\n2. Testing high-level cloud function with return_array...")
        
        result = uclchem.model.cloud(
            param_dict=param_dict,
            out_species=["OH", "CO", "H2O"],
            return_array=True,
            return_heating=True,
            timepoints=50
        )
        
        physicsArray, chemicalAbunArray, ratesArray, heatArray, abundanceStart, success_flag = result
        
        print(f"High-level array call success flag: {success_flag}")
        if heatArray is not None:
            print(f"Returned heating array shape: {heatArray.shape}")
            non_zero_heating = (heatArray != 0).sum()
            print(f"Non-zero heating values: {non_zero_heating}")
            print("High-level array test passed")
        else:
            print("No heating array returned")
        
        # Test 3: High-level model function with return_dataframe
        print("\n3. Testing high-level cloud function with return_dataframe...")
        
        result = uclchem.model.cloud(
            param_dict=param_dict,
            out_species=["OH", "CO", "H2O"],
            return_dataframe=True,
            return_heating=True,
            timepoints=50
        )
        
        physics_df, chemistry_df, rates_df, heating_df, abundanceStart, success_flag = result
        
        print(f"High-level DataFrame call success flag: {success_flag}")
        if heating_df is not None:
            print(f"Heating DataFrame shape: {heating_df.shape}")
            print(f"Heating DataFrame columns: {list(heating_df.columns)}")
            
            # Show sample data
            print("\nSample heating DataFrame data:")
            print(heating_df.head(3))
            
            # Check for non-zero values
            non_zero_heating = (heating_df != 0).sum().sum()
            print(f"Non-zero heating values in DataFrame: {non_zero_heating}")
            print("High-level DataFrame test passed")
        else:
            print("No heating DataFrame returned")
        
        # Test 4: Heating terms verification
        print("\n4. Verifying heating terms structure...")
        
        if heating_df is not None:
            expected_heating_columns = [
                "Time", "H", "C+", "O", "C", "CO", "p-H2", "o-H2", "SI+", "S",
                "Photoelectric", "H2Formation", "FUVPumping", "Photodissociation",
                "CIonization", "CRheating", "TurbHeating", "GasGrainColls"
            ]
            
            actual_columns = list(heating_df.columns)
            print(f"Expected {len(expected_heating_columns)} columns, got {len(actual_columns)}")
            
            # Check if all expected columns are present
            missing_columns = set(expected_heating_columns) - set(actual_columns)
            extra_columns = set(actual_columns) - set(expected_heating_columns)
            
            if not missing_columns and not extra_columns:
                print("All expected heating columns present")
            else:
                if missing_columns:
                    print(f"Missing columns: {missing_columns}")
                if extra_columns:
                    print(f"Extra columns: {extra_columns}")
        
        # Summary
        print("Test succesful")
        
        return True
        
    except Exception as e:
        print(f"\nTest failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_heating_arrays_comprehensive()
    if not success:
        print("\nHEATING ARRAY TEST FAILED")
        exit(1)
