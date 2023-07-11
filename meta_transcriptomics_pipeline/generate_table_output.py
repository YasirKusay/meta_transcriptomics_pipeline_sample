htmlCode = """
<html>
    <head>
        <title>My New Page!</title>
        <style>
            body {
                display: flex;
                margin: 0;
                flex-direction: column;
                min-height: 100%; /* to ensure we fill the remaining space*/
                font-family: sans-serif;
            }

            .title {
                padding: 10px;
                border: 5px;
                width: 100%;
                border-bottom: 1px gray;
            }

            .main_body {
                width: 100%;
                position: relative;
                padding-left: 20px;
                padding-right: 20px;
                overflow: auto;
                white-space: nowrap;
                display:flex;
                flex-direction: column;
                align-items: center;
                left: -20px;
                flex-grow: 1; /* used to fill the remaining space in the body */
            }

            .filtering_options {
                padding: 30px;
                width: 100%;
                display: flex;
                flex-direction: row;
                justify-content: space-evenly;
                flex-wrap: nowrap;
            }

            .autocomplete {
                width: 250px;
                position: absolute;
                background: white;
                top: 34px;
                z-index: 2;
            }

            .select_species {
                height: 25px;
                flex: 1;
                border-radius: 8px;
            }

            .arrow_down {
                position: relative;
                content: "";
                display: inline-block;
                width: 12px;
                height: 12px;
                border-right: 0.2em solid black;
                border-top: 0.2em solid black;
                transform: rotate(135deg);
                left: -25px;
                top: 2px;
            }

            .search_recommendation {
                padding: 5px;
                cursor: pointer;
            }

            .search_recommendation:not(:last-child) {
                border-bottom: 1px solid black;
            }

            .search_recommendation:hover {
                background-color: grey;
                width: 100%;
            }

            .select_species_dropdown {
                display: flex;
                flex-direction: column;
            }

            .toggle_dropdown_button {
                display: flex;
                flex-direction: row;
                flex-wrap: no-wrap;
            }

            .toggle_dropdown_button:hover {
                cursor: pointer;
            }

            .toggle_dropdown {
                height: 25px;
                border: 2px solid grey;
                width: 150px;
                border-radius: 8px;
                color: grey;
                padding-left: 2px;
                display: flex;
                align-items: center;
                -webkit-box-sizing: border-box; 
                -moz-box-sizing: border-box;    
                box-sizing: border-box;    
            }

            .select_species_options {
                width: 300px;
                height: 200px;
                visibility: hidden;
                position: absolute;
                display: flex;
                justify-content: space-between;
                top: 55px;
                background-color: rgba(179, 179, 179);
                border-radius: 8px;
                z-index: 2;
                padding: 5px;
            }

            .filter_by_value_options {
                width: 300px;
                height: 200px;
                visibility: hidden;
                position: absolute;
                top: 55px;
                background-color: rgba(179, 179, 179);
                border-radius: 8px;
                z-index: 2;
                padding: 5px;
                display: flex;
                justify-content: space-around;
                cursor: default;
            }

            .select_filtering_options {
                height: 25px;
                border: 2px solid grey;
                width: 125px;
                border-radius: 8px;
                color: grey;
                padding-left: 2px;
            }

            .apply_filter_button {
                height: 25px;
                width: 125px;
                border-radius: 8px;
            }
            
            #apply_filter_value.apply_filter_button {
                bottom: 10px;
                position: absolute;
            }

            .text_input {
                height: 25px;
                border: 2px solid grey;
                width: 125px;
                border-radius: 8px;
                color: grey;
                padding-left: 2px;
            }

            .applied_filters {
                width: 80%;
                display: flex;
                flex-direction: row;
                white-space: nowrap;
                justify-content: space-between;
            }

            .applied_filter {
                height: 30px;
                padding-left: 5px;
                padding-right: 5px;
                margin-bottom: 20px;
                border-radius: 8px;
                background: rgba(128, 128, 128, 0.8);
                filter: brightness(150%);
                font-size: 14px;
                display: flex;
                flex-direction: row;
                justify-content: space-evenly;
                align-items: center;
                cursor: default;
            }

            .invalid_input {
                color: red;
                visibility: hidden;
                position: absolute;
                top: 25px;
            }

            .result_box {
                background: rgba(128, 128, 128, 0.3);
                width: 95%;
                /* padding-left: 10px; */
                height: 30px;
                border-top: 1px solid gray;
                border-left: 1px solid gray;
                border-right: 1px solid gray;
                filter: brightness(150%);
                display: flex;
                flex-direction: row;
                align-items: center;
            }

            .lineage {
                display: none;
            }

            .result_box:last-child {
                border-bottom: 1px solid gray;
            }

            .result_element {
                display: flex;
                flex-direction: row;
                justify-content: center;
                overflow: hidden;
                white-space: nowrap;
            }

            .species_name {
                line-height: 30px;
                flex: 3;
            }

            .score_element {
                border-left: 1px solid gray;
                line-height: 30px;
                padding-left: 10px;
                flex: 1;
            }

            .result_header {
                filter: brightness(110%);
                border-bottom: 1px solid;
            }
            
            .header_field {
                position: static;
            }

            .header_field:hover {
                background: rgba(128, 128, 128, 0.3);
                filter: brightness(90%);
                cursor: pointer;
            }

            .species_name:hover {
                background: rgba(128, 128, 128, 0.3);
                filter: brightness(30%);
                cursor: pointer;
            }
        </style>
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">
    </head>

    <body>
        <div class="title">
            <h1>Output</h1>
        </div>
        <div class="main_body">
            <div class="filtering_options">
                <div class="select_species_dropdown">
                    <div class="toggle_dropdown_button" id="toggle_species_dropdown">
                        <div class="toggle_dropdown" id="toggle_species_dropdown_text"> 
                            Categories
                        </div>
                        <i class="arrow_down" id="select_species_dropdown_arrow"></i>
                    </div>
                    <div class="select_species_options">
                        <div class="filtering_textbox">
                            <input type="text" placeholder="Search.." class="text_input" id="lineage_to_filter_by">
                        </div>
                        <p class="invalid_input" id="display_input_warning">Invalid input.</p>
                        <button class="apply_filter_button" id="apply_filter_lineage">Apply Filter</button>
                        <div class="autocomplete"></div>
                    </div>
                </div>
                <div class="filter_by_value_dropdown">
                    <div class="toggle_dropdown_button" id="toggle_value_dropdown">
                        <div class="toggle_dropdown" id="toggle_filter_dropdown_text"> 
                            Filter Values
                        </div>
                        <i class="arrow_down" id="toggle_value_dropdown_arrow"></i>
                        <div class="filter_by_value_options">
                            <select class="select_filtering_options">
                                <option value="" disabled selected>Select your Filtering Option</option>
                                <option value="percentage_reads_mapped">Reads Mapped (%)</option>
                                <option value="tpm_score">TPM Score</option>
                                <option value="num_reads_mapped">Number of Reads Mapped</option>
                                <option value="average_e_value">Average E-value</option>
                                <option value="average_bitscore">Average Bitscore</option>
                                <option value="average_percentage_identity">Average Percentage Identity</option>
                                <option value="average_alignment_length">Average Alignment Length</option>
                            </select>
                            <div class="filtering_textbox">
                                <input type="text" placeholder="Search.." class="text_input" id="value_to_filter_by">
                            </div>
                            <p class="invalid_input" id="display_input_warning">Invalid input.</p>
                            <button id="apply_filter_value" class="apply_filter_button">Apply Filter</button>
                        </div>
                    </div>
                </div>
            </div>
            <div class="applied_filters">

            </div>
            <div class="result_box result_header">
                <div class="result_element species_name header_field" id="organism_header">
                    Organism
                </div>
                <div class="result_element score_element header_field" id="percentage_reads_mapped">
                    Reads mapped (%)
                </div>
                <div class="result_element score_element header_field" id="tpm_score">
                    TPM score
                </div>
                <div class="result_element score_element header_field" id="num_reads_mapped">
                    Reads mapped
                </div>
                <div class="result_element score_element header_field" id="average_e_value">
                    E-value
                </div>
                <div class="result_element score_element header_field" id="average_bitscore">
                    Bitscore
                </div>
                <div class="result_element score_element header_field" id="average_percentage_identity">
                    Percentage Id
                </div>
                <div class="result_element score_element header_field" id="average_alignment_length">
                    Avg Alignment Length
                </div>
            </div>
"""

jsScript = """
            var all_lineages = [] // stores all unique lineages
            // adding all unique lineages to all_lineages:
            var all_displayed_species = document.getElementsByClassName("result_box")
            for (let i = 0; i < all_displayed_species.length; i++) {
                var curr_species = all_displayed_species[i]
                if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;
                var curr_lineage = curr_species.getElementsByClassName('lineage')[0].innerHTML.trim().split(",")
                for (let j = 0; j < curr_lineage.length; j++) {
                    if (curr_lineage[j].trim().toLowerCase() == "unknown") continue;
                    if (!all_lineages.includes(curr_lineage[j].trim().toLowerCase())) all_lineages.push(curr_lineage[j].trim().toLowerCase())
                }
            }

            var lineages_filtered_list = []

            // determining the position for the e_value within the header
            var e_value_position = -1
            var filtering_value_headers = document.getElementsByClassName('result_header')[0].getElementsByTagName("div");
            for (var i = 1; i < filtering_value_headers.length; i++) {
                if (filtering_value_headers[i].id === "average_e_value") {
                    e_value_position = i;
                    break;
                }
            }


            var select_species_dropdown_elem = document.getElementsByClassName('select_species_dropdown')[0]
            var filter_by_value_dropdown_elem = document.getElementsByClassName('filter_by_value_dropdown')[0]

            var select_species_options_elem = document.getElementsByClassName('select_species_options')[0]
            var filter_by_value_options_elem = document.getElementsByClassName('filter_by_value_options')[0]

            var all_species_name = document.getElementsByClassName('species_name')
            var autocomplete_element = document.getElementsByClassName('autocomplete')[0]
            var invalid_input_elem = document.getElementById('display_input_warning')

            select_species_dropdown_elem.addEventListener('click', () => {select_species_options_elem.style.visibility = 'visible';})
            filter_by_value_dropdown_elem.addEventListener('click', () => {filter_by_value_options_elem.style.visibility = 'visible';})

            // hashmap of lists
            // lust contains things currently filtered out for each filtering option
            var currently_filtered_species = new Map([["lineage", []],
                                                    ["percentage_reads_mapped", []], 
                                                    ["tpm_score", []], 
                                                    ["num_reads_mapped", []], 
                                                    ["average_e_value", []], 
                                                    ["average_bitscore", []], 
                                                    ["average_percentage_identity", []],
                                                    ["average_alignment_length", []]]
                                                    );

            // skip refers to what column we dont want to check
            // typically will be set to be equal to the category 
            // that is being filtered/unfiltered
            function check_if_currently_filtered(species, skip) {
            for (var currFilteringOption of currently_filtered_species.keys()) {
                if (currFilteringOption === skip) continue;
                if (currently_filtered_species.get(currFilteringOption).includes(species)) return true;
            }
            return false;
            }

            window.addEventListener('click', ({ target }) => {
            if (!select_species_dropdown_elem.contains(target)) select_species_options_elem.style.visibility = 'hidden';
            if (!filter_by_value_dropdown_elem.contains(target)) {
                filter_by_value_options_elem.style.visibility = 'hidden';
                invalid_input_elem.style.visibility = "hidden";
            }
            });

            var lineage_to_filter_by = document.getElementById('lineage_to_filter_by');
            lineage_to_filter_by.addEventListener('keyup', suggestSpecies)

            function capitalise(string) {
                return string[0].toUpperCase() + string.slice(1);
            }

            function suggestSpecies() {
                var searched = lineage_to_filter_by.value.trim().toLowerCase();
                autocomplete_element.innerHTML = ""
                if (searched.length < 4) return;

                var toAdd = []

                // need to grab the species
                // skipping first element, as it is the title
                for (let i = 0; i < all_lineages.length; i++) {
                    let curr_lineage = all_lineages[i].trim();
                    if (curr_lineage.includes(searched) === true) toAdd.push(curr_lineage);
                    if (toAdd.length === 10) break // do not display more than 10 elements
                }

                for (let i = 0; i < toAdd.length; i++) {
                    var search_result_elem = document.createElement('div');
                    search_result_elem.className = "search_recommendation"
                    search_result_elem.innerHTML = toAdd[i].trim().replace('\\n','')
                    autocomplete_element.appendChild(search_result_elem)
                    search_result_elem.addEventListener('click', ({target}) => {
                        lineage_to_filter_by.value = target.innerHTML;
                        autocomplete_element.innerHTML = ""
                        add_filter_by_lineage();
                    })
                } 
            }

            var filter_area = document.getElementsByClassName('applied_filters')[0]

            function clear_lineage_filtering({ target }) {
                var filter_elem = (target.parentElement.className === "applied_filter") ? target.parentElement : target;
                // need to remove the lineage from lineages_filtered_list
                var lineage_to_remove = filter_elem.id.split("_")[0].trim().toLowerCase() // format is (lineage)_filter_active,
                remove_elem_from_list(lineages_filtered_list, lineage_to_remove)
                filter_elem.remove();
                check_filter_by_lineage();
            }

            // need a way to ignore case
            function add_filter_by_lineage() {
                var selected_lineage = lineage_to_filter_by.value.trim().toLowerCase();
                if (lineages_filtered_list.includes(selected_lineage)) return; // exit if lineage in currently_active_lineages
                if (!all_lineages.includes(selected_lineage)) return;

                // we need to check if our lineage is a valid lineage
                // we should also take the case of something into account

                lineages_filtered_list.push(selected_lineage) // add lineage to currently_active_lineages

                // add a visual filtering "bubble"
                var filter_display = document.createElement('div');
                filter_display.className = "applied_filter";
                var text = document.createElement('p');
                text.innerHTML = capitalise(selected_lineage);
                filter_display.id = selected_lineage.concat("_filter_active");
                filter_display.appendChild(text);
                filter_display.addEventListener('click', clear_lineage_filtering)
                filter_area.appendChild(filter_display)

                check_filter_by_lineage();
            }

            document.getElementById('apply_filter_lineage').addEventListener('click', add_filter_by_lineage)

            function check_filter_by_lineage() {
                var all_displayed_species = document.getElementsByClassName("result_box")
                // nothing selected, remove every species from currently_filtered_species['lineage']
                if (lineages_filtered_list.length == 0) {
                    for (let i = 0; i < all_displayed_species.length; i++) {
                        let curr_species = all_displayed_species[i]
                        if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;
                        // for this scenario, its better to use capitalisation in the currently_filtered_species, so we match filter by value
                        if (!check_if_currently_filtered(curr_species.getElementsByClassName('species_name')[0].innerHTML.trim(), 'lineage')) curr_species.style.display = "";
                    }
                    currently_filtered_species.set('lineage', [])
                } else {
                    for (let i = 0; i < all_displayed_species.length; i++) {
                        let curr_species = all_displayed_species[i]
                        var curr_species_name = curr_species.getElementsByClassName('species_name')[0].innerHTML.trim()
                        if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;

                        // need to determine if any lineages in lineages_filtered_list are part of the lineages of our current species
                        // i.e. if there is an intersection between lineages_filtered_list and the lineages of the current species
                        var curr_lineage = curr_species.getElementsByClassName('lineage')[0].innerHTML.trim().split(",")
                        for (let j = 0; j < curr_lineage.length; j++) {
                            curr_lineage[j] = curr_lineage[j].trim().toLowerCase() // to lowercase to match lineages_filtered_list
                        }

                        var intersection = lineages_filtered_list.filter(value => curr_lineage.includes(value));
                        if (intersection.length > 0) {
                            if (!check_if_currently_filtered(curr_species_name, 'lineage')) curr_species.style.display = "";
                            currently_filtered_species.set('lineage', remove_elem_from_list(currently_filtered_species.get('lineage'), curr_species_name))
                        } else {
                            var lineage_filtered_list = currently_filtered_species.get("lineage");
                            if (!lineage_filtered_list.includes(curr_species_name)) {
                                lineage_filtered_list.push(curr_species_name);
                                currently_filtered_species.set('lineage', lineage_filtered_list)
                            }
                            curr_species.style.display = "none";
                        }
                    }
                }
            }

            var submit_value_filter = document.getElementById('apply_filter_value')
            var select_filtering_options_elem = document.getElementsByClassName('select_filtering_options')[0]
            var value_to_filter_by_elem = document.getElementById('value_to_filter_by')

            function updateFilterElemIfActive(elem, name, comparison_sign, value) {
            if (elem === null) return false;
            elem.innerHTML = name + " " + comparison_sign + " " + value;
            return true;
            }

            // current method is very flawed, what if order of select_filtering_options doesnt match our header order
            function get_index_of_select_element(name) {
            var filtering_value_options = document.getElementsByClassName('result_header')[0].getElementsByTagName("div");
            var index_of_selected_elem = -1;
            for (var i= 1; i < filtering_value_options.length; i++) {
                if (filtering_value_options[i].id === name) {
                    index_of_selected_elem = i;
                    break;
                }
            }
            return index_of_selected_elem;
            }

            function clear_value_filter({ target }) {
                var filter_elem = (target.parentElement.className === "applied_filter") ? target.parentElement : target;
                // we can get the name of the select tag from target.id (just need to trim _active from it)
                // using that, we can get the index of the select elements.
                var name_of_select_elem = filter_elem.id.replace('_active', '');
                // var value_filtered_by

                var all_displayed_species = document.getElementsByClassName("result_box")
                for (let i = 0; i < all_displayed_species.length; i++) {
                    let curr_species = all_displayed_species[i]
                    var curr_species_name = curr_species.getElementsByClassName('species_name')[0].innerHTML.trim()
                    if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;
                    if (!check_if_currently_filtered(curr_species_name, name_of_select_elem)) curr_species.style.display = "";
                }
                currently_filtered_species.set(name_of_select_elem, [])

                console.log(currently_filtered_species)
                // now we can reverse the select 
                filter_elem.remove();
            } 

            function click_apply_value_filter() {
            invalid_input_elem.style.visibility = "hidden";
            if (
                (select_filtering_options_elem.value === "average_e_value" && value_to_filter_by_elem.value > 1) ||
                ((select_filtering_options_elem.value === "percentage_reads_mapped" || select_filtering_options_elem.value === "average_percentage_identity") && value_to_filter_by_elem.value >= 100) ||
                (value_to_filter_by_elem.value < 0)
                ) {
                value_to_filter_by_elem.value = "";
                invalid_input_elem.style.visibility = "visible";
                // also throw a suitable warning
                return;
            }

            // now lets apply the popup here
            // but if it already exists, then update it
            var comparison_sign = ">";
            var filter_value_name = "";
            var filter_value_name_id = "";
            if (select_filtering_options_elem.value === "percentage_reads_mapped") {
                filter_value_name = "Reads Mapped (%)";
                filter_value_name_id = "percentage_reads_mapped_active";
            }
            if (select_filtering_options_elem.value === "tpm_score") {
                filter_value_name = "TPM Score";
                filter_value_name_id = "tpm_score_active";
            }
            if (select_filtering_options_elem.value === "num_reads_mapped") {
                filter_value_name = "Num Reads Mapped";
                filter_value_name_id = "num_reads_mapped_active";
            }
            if (select_filtering_options_elem.value === "average_e_value") {
                filter_value_name = "Average E-value";
                comparison_sign = "<";
                filter_value_name_id = "average_e_value_active";
            }
            if (select_filtering_options_elem.value === "average_bitscore") {
                filter_value_name = "Average Bitscore";
                filter_value_name_id = "average_bitscore_active";
            }
            if (select_filtering_options_elem.value === "average_percentage_identity") {
                filter_value_name = "Average Percentage Identity";
                filter_value_name_id = "average_percentage_identity_active";
            }
            if (select_filtering_options_elem.value === "average_alignment_length") {
                filter_value_name = "Average Alignment Length";
                filter_value_name_id = "average_alignment_length_active";
            }

            // if (updateFilterElemIfActive(document.getElementById(filter_value_name_id), filter_value_name, comparison_sign, value_to_filter_by_elem.value)) return;
            if (!updateFilterElemIfActive(document.getElementById(filter_value_name_id), filter_value_name, comparison_sign, value_to_filter_by_elem.value)) {
                var filter_display = document.createElement('div');
                filter_display.className = "applied_filter";
                var text = document.createElement('p');
                text.innerHTML = filter_value_name + " " + comparison_sign + " " + parseFloat(value_to_filter_by_elem.value);
                filter_display.id = filter_value_name_id;
                filter_display.style.width = "200px";
                filter_display.appendChild(text)
                filter_display.addEventListener('click', clear_value_filter)
                filter_area.appendChild(filter_display)
            }

            // now we need to actually filter the tables.
            index_of_selected_elem = get_index_of_select_element(select_filtering_options_elem.value)

            // the above calculation includes the select your filtering option "choice", need to
            // take that into account

            // now we need to grab the result_box classes, we can grab the corresponding element we want to compare 
            // don't forget to remove the header tag
            // also the index of E-value cell is 4, need to see if the value of e-value of cell is less than our
            // entered e-value, and to keep it if it is
            var all_displayed_species = document.getElementsByClassName("result_box")
            for (let i = 0; i < all_displayed_species.length; i++) {
                let curr_species = all_displayed_species[i]
                // curr_species.style.display = ""; bad idea, will revert any previously applied filters
                if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;
                // first comparison is the e-value case

                //if ((index_of_selected_elem === 4 && (parseFloat(curr_species.getElementsByTagName('div')[index_of_selected_elem].innerHTML.trim()) > parseFloat(value_to_filter_by_elem.value.trim()))) ||
                //(index_of_selected_elem !== 4 && (parseFloat(curr_species.getElementsByTagName('div')[index_of_selected_elem].innerHTML.trim()) < parseFloat(value_to_filter_by_elem.value.trim())))) {
                
                // below is actually a better way of implementing which column was selected, makes it much easier to include new headers
                if ((select_filtering_options_elem.value === "average_e_value" && (parseFloat(curr_species.getElementsByTagName('div')[index_of_selected_elem].innerHTML.trim()) > parseFloat(value_to_filter_by_elem.value.trim()))) ||
                (select_filtering_options_elem.value !== "average_e_value" && (parseFloat(curr_species.getElementsByTagName('div')[index_of_selected_elem].innerHTML.trim()) < parseFloat(value_to_filter_by_elem.value.trim())))) {
                    curr_species.style.display = "none";
                    var value_filtered_list = currently_filtered_species.get(select_filtering_options_elem.value);
                    var curr_species_name = curr_species.getElementsByClassName('species_name')[0].innerHTML.trim()
                    if (!value_filtered_list.includes(curr_species_name)) {
                    value_filtered_list.push(curr_species_name);
                    currently_filtered_species.set(select_filtering_options_elem.value, value_filtered_list)
                    }
                    curr_species.style.display = "none";
                }
            }

            // now need to work out removing filters
            // when removing the elements, we need to still check all currently active filters to determine if we 
            // can revert an entry

            console.log(currently_filtered_species)
            }

            submit_value_filter.addEventListener('click', click_apply_value_filter);

            function remove_elem_from_list(list, value) {
            for (let i = 0; i < list.length; i++) {
                if (list[i].trim() === value.trim()) {
                list.splice(i, 1);
                break;
                }
            }
            return list;
            }

            function sortFloat(a,b) { return a - b; }

            organism_header = document.getElementById("organism_header")
            organism_header.addEventListener('click', () => {reorderTable(0)})

            perct_reads_mapped_header = document.getElementById("percentage_reads_mapped")
            perct_reads_mapped_header.addEventListener('click', () => {reorderTable(get_index_of_select_element("percentage_reads_mapped"))})

            tpm_score_header = document.getElementById("tpm_score")
            tpm_score_header.addEventListener('click', () => {reorderTable(get_index_of_select_element("tpm_score"))})

            num_reads_mapped_header = document.getElementById("num_reads_mapped")
            num_reads_mapped_header.addEventListener('click', () => {reorderTable(get_index_of_select_element("num_reads_mapped"))})

            e_val_header = document.getElementById("average_e_value")
            e_val_header.addEventListener('click', () => {reorderTable(get_index_of_select_element("average_e_value"))})

            bitscore_header = document.getElementById("average_bitscore")
            bitscore_header.addEventListener('click', () => {reorderTable(get_index_of_select_element("average_bitscore"))})

            perct_identity_header = document.getElementById("average_percentage_identity")
            perct_identity_header.addEventListener('click', () => {reorderTable(get_index_of_select_element("average_percentage_identity"))})

            avg_length_header = document.getElementById("average_alignment_length")
            avg_length_header.addEventListener('click', () => {reorderTable(get_index_of_select_element("average_alignment_length"))})

            // used to store the previously selected column to sort by
            // if what you select now is equal to this, sort it in reverse
            current_sort_selected = -1

            function reorderTable(position_of_elem) {
                if (current_sort_selected !== position_of_elem) {
                    if (position_of_elem !== 0) {
                        var species_value_map = new Map();
                        var all_displayed_species = document.getElementsByClassName("result_box")
                        var sorted = []
                        for (let i = 0; i < all_displayed_species.length; i++) {
                            let curr_species = all_displayed_species[i]
                            if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;
                            species_value_map.set(curr_species.getElementsByTagName('div')[0].innerHTML.trim(), // name of element
                                parseFloat(curr_species.getElementsByTagName('div')[position_of_elem].innerHTML.trim())) // value of element
                        }

                        // sorting hashmap by key 
                        // from https://stackoverflow.com/questions/37982476/how-to-sort-a-map-by-value-in-javascript
                        
                        // e_value will be ordered by descending order
                        var map_sorted;
                        if (position_of_elem !== e_value_position)  {
                            map_sorted = new Map([...species_value_map.entries()].sort((a, b) => b[1] - a[1]));
                        } else {
                            map_sorted = new Map([...species_value_map.entries()].sort((a, b) => a[1] - b[1]));
                        }

                        const sorted_names = Array.from(map_sorted.keys())

                        for (let i = 0; i < all_displayed_species.length; i++) {
                            let curr_species = all_displayed_species[i]
                            if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;
                            console.log(sorted_names.indexOf(curr_species.getElementsByTagName('div')[0].innerHTML.trim()))
                            curr_species.style.order = sorted_names.indexOf(curr_species.getElementsByTagName('div')[0].innerHTML.trim()) + 1
                            // change the thing at the bottom with the black border also
                            
                            // last species
                            if (curr_species.style.order == all_displayed_species.length - 1) {
                                curr_species.style.borderBottom = "1px solid gray";
                            } else {
                                curr_species.style.borderBottom = "none";
                            }
                        }
                    } else { // we are sorting the names of the organisms alphabetically
                        const sorted_names = []

                        var all_displayed_species = document.getElementsByClassName("result_box")
                        for (let i = 0; i < all_displayed_species.length; i++) {
                            let curr_species = all_displayed_species[i]
                            if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;
                            sorted_names.push(curr_species.getElementsByTagName('div')[0].innerHTML.trim());
                        }

                        sorted_names.sort();

                        for (let i = 0; i < all_displayed_species.length; i++) {
                            console.log(all_displayed_species.length)
                            let curr_species = all_displayed_species[i]
                            if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;
                            console.log(curr_species.getElementsByTagName('div')[0].innerHTML.trim())
                            console.log(sorted_names.indexOf(curr_species.getElementsByTagName('div')[0].innerHTML.trim()))
                            curr_species.style.order = sorted_names.indexOf(curr_species.getElementsByTagName('div')[0].innerHTML.trim()) + 1
                            // change the thing at the bottom with the black border also
                            
                            // last species
                            if (curr_species.style.order == all_displayed_species.length - 1) {
                                curr_species.style.borderBottom = "1px solid gray";
                            } else {
                                curr_species.style.borderBottom = "none";
                            }
                        }
                    }
                } else {
                    console.log("REVERSE")
                    // sorting the previously selected column but in reverse
                    // simply reverse what we have right now
                    var all_displayed_species = document.getElementsByClassName("result_box")
                    for (let i = 0; i < all_displayed_species.length; i++) {
                        let curr_species = all_displayed_species[i]
                        if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;

                        let old = curr_species.style.order
                        // below may not make sense but it is important for maintaining the order of the elements
                        // if it was all_displayed_species - 1 - order, the order could become 0
                        curr_species.style.order = all_displayed_species.length - curr_species.style.order
                        console.log(old, curr_species.style.order)
                        // change the thing at the bottom with the black border also
                        
                        // last species
                        if (curr_species.style.order == all_displayed_species.length - 1) {
                            curr_species.style.borderBottom = "1px solid gray";
                        } else {
                            curr_species.style.borderBottom = "none";
                        }
                    }
                }
                current_sort_selected = position_of_elem
            }
"""

def generate_table_output(tableSummaryInput, outputHTMLSummary):
    global htmlCode
    with open(tableSummaryInput, "r") as f:
        for line in f:
            curr = line.split("\t") 
            organism = curr[0]
            percent_reads_mapped = curr[1]
            tpm_score = curr[2]
            num_reads_mapped = curr[3]
            e_value = curr[4]
            bitscore = curr[5]
            avg_percent_identity = curr[6]
            avg_alignment_len = curr[7]
            lineage = curr[8].strip() # final column
            htmlCode += "            <div class=\"result_box\">\n"
            htmlCode += "                <div class=\"result_element species_name\">\n"
            htmlCode += "                    " + str(organism) + "\n"
            htmlCode += "                </div>\n"

            for score in [percent_reads_mapped, tpm_score, num_reads_mapped, e_value, bitscore, avg_percent_identity, avg_alignment_len]:
                htmlCode += "                <div class=\"result_element score_element\">\n"
                htmlCode += "                    " + str(score) + "\n"
                htmlCode += "                </div>\n"

            htmlCode += "                <div class=\"lineage\">\n"
            htmlCode += "                    " + lineage + "\n"
            htmlCode += "                </div>\n"

            htmlCode += "            </div>\n"

    htmlCode += "        </div>\n"
    htmlCode += "        <script>\n"
    htmlCode += jsScript
    htmlCode += "        </script>\n"
    htmlCode += "    </body>\n"
    htmlCode += "</html>\n"
    
    wf = open(outputHTMLSummary, "w")
    wf.write(htmlCode)
    wf.close()
