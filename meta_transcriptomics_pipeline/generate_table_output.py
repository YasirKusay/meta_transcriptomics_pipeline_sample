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

            .select_species {
                height: 25px;
                flex: 1;
                border-radius: 8px;
            }

            .search_species {
                display: flex;
                flex-direction: column;
                position: relative;
            }

            .search_species_input {
                height: 25px;
                width: 150px;
                border-radius: 8px;
                border-color: grey;
            }

            .autocomplete {
                width: 150px;
                position: absolute;
                background: gray;
                filter: brightness(150%);
                top: 25px;
                z-index: 2;
            }

            .submit_search_species {
                position: relative;
                left: -32px;
                background: white;
                border: none;
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
                width: 150px;
                visibility: hidden;
                position: absolute;
                top: 55px;
                background-color: rgba(128, 128, 128);
                filter: brightness(150%);
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

            .apply_filter {
                position: absolute;
                bottom: 10px;
                height: 25px;
                width: 125px;
                border-radius: 8px;
            }

            .value_to_filter_by {
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
                width: 100px;
                margin-bottom: 20px;
                border-radius: 8px;
                background: rgba(64, 224, 208, 1.0);
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
                padding-left: 10px;
                height: 30px;
                border-top: 1px solid gray;
                border-left: 1px solid gray;
                border-right: 1px solid gray;
                filter: brightness(150%);
                display: flex;
                flex-direction: row;
                align-items: center;
            }

            .domain {
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
            }

            .species_name {
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
        </style>
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">
    </head>

    <body>
        <div class="title">
            <h1>Output</h1>
        </div>
        <div class="main_body">
            <div class="filtering_options">
                <div class="search_species">
                    <div class="search_field">
                        <input type="text" placeholder="Search.." class="search_species_input">
                        <button type="submit" class="submit_search_species"><i class="fa fa-search search_button"></i>
                    </div>
                    <div class="autocomplete"></div>
                </div>
                <div class="select_species_dropdown">
                    <div class="toggle_dropdown_button" id="toggle_species_dropdown">
                        <div class="toggle_dropdown" id="toggle_species_dropdown_text"> 
                            Categories
                        </div>
                        <i class="arrow_down" id="select_species_dropdown_arrow"></i>
                    </div>
                    <div class="select_species_options">
                        <div class="filter_species">
                            <input type="checkbox" id="eukaryotes_filter">
                            Eukaryotes
                        </div>
                        <div class="filter_species">
                            <input type="checkbox" id="bacteria_filter">
                            Bacteria
                        </div>
                        <div class="filter_species">
                            <input type="checkbox" id="virus_filter">
                            Viruses
                        </div>
                        <div class="filter_species">
                            <input type="checkbox" id="archaea_filter">
                            Archaea
                        </div>
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
                                <option value="average_query_length">Average Query Length</option>
                            </select>
                            <div class="filtering_textbox">
                                <input type="text" placeholder="Search.." class="value_to_filter_by">
                            </div>
                            <p class="invalid_input" id="display_input_warning">Invalid input.</p>
                            <button class="apply_filter">Apply Filter</button>
                        </div>
                    </div>
                </div>
            </div>
            <div class="applied_filters">

            </div>
            <div class="result_box result_header">
                <div class="result_element species_name">
                    Organism
                </div>
                <div class="result_element score_element">
                    Reads mapped
                </div>
                <div class="result_element score_element">
                    TPM score
                </div>
                <div class="result_element score_element">
                    Num reads mapped
                </div>
                <div class="result_element score_element">
                    E-value
                </div>
                <div class="result_element score_element">
                    Bitscore
                </div>
                <div class="result_element score_element">
                Percentage Identity (%)
                </div>
                <div class="result_element score_element">
                    Average Length
                </div>
            </div>
"""

jsScript = """
            var select_species_dropdown_elem = document.getElementsByClassName('select_species_dropdown')[0]
            var filter_by_value_dropdown_elem = document.getElementsByClassName('filter_by_value_dropdown')[0]

            var select_species_options_elem = document.getElementsByClassName('select_species_options')[0]
            var filter_by_value_options_elem = document.getElementsByClassName('filter_by_value_options')[0]

            var search_species_input_dropdown = document.getElementsByClassName('search_species_input')[0]
            var all_species_name = document.getElementsByClassName('species_name')
            var autocomplete_element = document.getElementsByClassName('autocomplete')[0]
            var invalid_input_elem = document.getElementById('display_input_warning')

            select_species_dropdown_elem.addEventListener('click', () => {select_species_options_elem.style.visibility = 'visible';})
            filter_by_value_dropdown_elem.addEventListener('click', () => {filter_by_value_options_elem.style.visibility = 'visible';})

            search_species_input_dropdown.addEventListener('keyup', suggestSpecies)

            // hashmap of lists
            // lust contains things currently filtered out for each filtering option
            var currently_filtered_species = new Map([["domain", []],
                                                    ["percentage_reads_mapped", []], 
                                                    ["tpm_score", []], 
                                                    ["num_reads_mapped", []], 
                                                    ["average_e_value", []], 
                                                    ["average_bitscore", []], 
                                                    ["average_percentage_identity", []],
                                                    ["average_query_length", []]]
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

            function suggestSpecies({}) {
            var searched = search_species_input_dropdown.value;

            autocomplete_element.innerHTML = ""

            if (searched.length < 4) {
                return;
            }

            var toAdd = []

            // need to grab the species
            // skipping first element, as it is the title
            for (let i = 1; i < all_species_name.length; i++) {
                let curr_species = all_species_name[i].innerHTML;
                if (curr_species.includes(searched) === true) {
                toAdd.push(curr_species)
                }
            }

            var toConcat = ""
            for (let i = 0; i < toAdd.length; i++) {
                toConcat = toConcat + toAdd[i].trim().replace('\\n','') + "</br>\\n"
            } 

            console.log(toConcat)

            autocomplete_element.innerHTML = toConcat;
            }


            window.addEventListener('click', ({ target }) => {
            if (!select_species_dropdown_elem.contains(target)) select_species_options_elem.style.visibility = 'hidden';
            if (!filter_by_value_dropdown_elem.contains(target)) {
                filter_by_value_options_elem.style.visibility = 'hidden';
                invalid_input_elem.style.visibility = "hidden";
            }
            });


            var filter_area = document.getElementsByClassName('applied_filters')[0]

            var eukaryotes_check_box = document.getElementById('eukaryotes_filter')
            var bacteria_filter_box = document.getElementById('bacteria_filter')
            var virus_filter_box = document.getElementById('virus_filter')
            var archaea_filter_box = document.getElementById('archaea_filter')

            function clear_filter_by_name(elem) {
            if (elem.id === "eukaryotes_filter_active") var temp = document.getElementById('eukaryotes_filter').checked = false
            if (elem.id === "bacteria_filter_active") var temp = document.getElementById('bacteria_filter').checked = false
            if (elem.id === "virus_filter_active") var temp = document.getElementById('virus_filter').checked = false
            if (elem.id === "archaea_filter_active") var temp = document.getElementById('archaea_filter').checked = false
            }

            function clear_filtering_checkbox({ target }) {
            if (target.parentElement.className === "applied_filter") {
                // also need to uncheck it from the checkbox
                // get the id to determine the correct filter
                clear_filter_by_name(target.parentElement)
                target.parentElement.remove()
            } else if (target.className === "applied_filter") {
                clear_filter_by_name(target)
                target.remove()
            }
            filter_by_domain()
            } 

            function species_filter_checkbox({ target }) {
            if (target.checked) {
                // add a "visual" filtering tag 
                var filter_display = document.createElement('div');
                filter_display.className = "applied_filter";
                var text = document.createElement('p');
                if (target.id === "eukaryotes_filter") {
                text.innerHTML = "Eukaryotes";
                filter_display.id = "eukaryotes_filter_active";
                }
                if (target.id === "bacteria_filter") {
                text.innerHTML = "Bacteria";
                filter_display.id = "bacteria_filter_active";
                }
                if (target.id === "virus_filter") { 
                text.innerHTML = "Viruses";
                filter_display.id = "virus_filter_active";
                }
                if (target.id === "archaea_filter") { 
                text.innerHTML = "Archaea";
                filter_display.id = "archaea_filter_active";
                }
                var x_icon = document.createElement('p');
                filter_display.appendChild(text)
                filter_display.addEventListener('click', clear_filtering_checkbox)
                filter_area.appendChild(filter_display)

                // now need to update the table according to the active species filter
            } else {
                if (target.id === "eukaryotes_filter") {
                var toRemove = document.getElementById("eukaryotes_filter_active").remove()
                console.log("hehe")
                }
                if (target.id === "bacteria_filter") {
                var toRemove = document.getElementById("bacteria_filter_active").remove()
                }
                if (target.id === "virus_filter") {
                var toRemove = document.getElementById("virus_filter_active").remove()
                }
                if (target.id === "archaea_filter") {
                var toRemove = document.getElementById("archaea_filter_active").remove()
                }
                // remove the visual filtering tag
            }
            filter_by_domain()
            }

            eukaryotes_check_box.addEventListener('click', species_filter_checkbox)
            bacteria_filter_box.addEventListener('click', species_filter_checkbox)
            virus_filter_box.addEventListener('click', species_filter_checkbox)
            archaea_filter_box.addEventListener('click', species_filter_checkbox)

            var submit_filter = document.getElementsByClassName('apply_filter')[0]
            var select_filtering_options_elem = document.getElementsByClassName('select_filtering_options')[0]
            var value_to_filter_by_elem = document.getElementsByClassName('value_to_filter_by')[0]

            function updateFilterElemIfActive(elem, name, comparison_sign, value) {
            if (elem === null) return false;
            elem.innerHTML = name + " " + comparison_sign + " " + value;
            return true;
            }

            function get_index_of_select_element(name) {
            var filtering_value_options = document.getElementsByClassName('select_filtering_options')[0].options;
            var index_of_selected_elem = -1;
            for (var i= 0; i < filtering_value_options.length; i++) {
                if (filtering_value_options[i].value === name) {
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

            function click_apply_filter() {
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
            console.log(document.getElementById("percentage_reads_mapped_active"));
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
            if (select_filtering_options_elem.value === "average_query_length") {
                filter_value_name = "Average Query Length";
                filter_value_name_id = "average_query_length_active";
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
                if ((index_of_selected_elem === 4 && (parseFloat(curr_species.getElementsByTagName('div')[index_of_selected_elem].innerHTML.trim()) > parseFloat(value_to_filter_by_elem.value.trim()))) ||
                (index_of_selected_elem !== 4 && (parseFloat(curr_species.getElementsByTagName('div')[index_of_selected_elem].innerHTML.trim()) < parseFloat(value_to_filter_by_elem.value.trim())))) {
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

            submit_filter.addEventListener('click', click_apply_filter);

            function remove_elem_from_list(list, value) {
            for (let i = 0; i < list.length; i++) {
                if (list[i].trim() === value.trim()) {
                list.splice(i, 1);
                break;
                }
            }
            return list;
            }

            function filter_by_domain() {
            var active_domain_filters = [];
            if (document.getElementById('eukaryotes_filter').checked === true) active_domain_filters.push("Eukaryotes");
            if (document.getElementById('bacteria_filter').checked === true) active_domain_filters.push("Bacteria");
            if (document.getElementById('virus_filter').checked === true) active_domain_filters.push("Viruses");
            if (document.getElementById('archaea_filter').checked === true) active_domain_filters.push("Archaea");

            // ensure everything is visible (unless other filters are applied on them)
            // IF THE ELEMENT IS CURRENTLY FILTERED OUT, PUT IT IN currently_filtered_species

            var all_displayed_species = document.getElementsByClassName("result_box")
            if (active_domain_filters.length === 0) {
                for (let i = 0; i < all_displayed_species.length; i++) {
                let curr_species = all_displayed_species[i]
                if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;
                if (!check_if_currently_filtered(curr_species.getElementsByClassName('species_name')[0].innerHTML.trim(), 'domain')) curr_species.style.display = "";
                currently_filtered_species.set('domain', [])
                }
            } else {
                for (let i = 0; i < all_displayed_species.length; i++) {
                curr_species = all_displayed_species[i]
                var curr_species_name = curr_species.getElementsByClassName('species_name')[0].innerHTML.trim()
                if (curr_species.contains(document.getElementsByClassName("result_header")[0])) continue;

                if (active_domain_filters.includes(curr_species.getElementsByClassName('domain')[0].innerHTML.trim())) {
                    if (!check_if_currently_filtered(curr_species_name, 'domain')) curr_species.style.display = "";
                    // remove it from here anyway
                    // this process will get quite confusing
                    currently_filtered_species.set('domain', remove_elem_from_list(currently_filtered_species.get('domain'), curr_species_name))
                } else {
                    var domain_filtered_list = currently_filtered_species.get("domain");
                    if (!domain_filtered_list.includes(curr_species_name)) {
                    domain_filtered_list.push(curr_species_name);
                    currently_filtered_species.set('domain', domain_filtered_list)
                    }
                    curr_species.style.display = "none";
                }
                }
            }
            console.log(currently_filtered_species)
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
            avg_query_len = curr[7]
            domain = curr[8].strip() # final column
            htmlCode += "            <div class=\"result_box\">\n"
            htmlCode += "                <div class=\"result_element species_name\">\n"
            htmlCode += "                    " + str(organism) + "\n"
            htmlCode += "                </div>\n"

            for score in [percent_reads_mapped, tpm_score, num_reads_mapped, e_value, bitscore, avg_percent_identity, avg_query_len]:
                htmlCode += "                <div class=\"result_element score_element\">\n"
                htmlCode += "                    " + str(score) + "\n"
                htmlCode += "                </div>\n"

            htmlCode += "                <div class=\"domain\">\n"
            htmlCode += "                    " + domain + "\n"
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
