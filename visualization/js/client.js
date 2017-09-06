//client.js : JavaScript code to that contacts server.

// Executes HMM stack and calls visualization with return aligned json
// stack: the stack to align as a json object
// rc: a boolean representing whether to use radiocarbon data
function run_hmm_stack(stack, rc) {

    // A new XMLHttpRequest object
    var request = new XMLHttpRequest();

    //Use MPS RESTful API to specify URL
    var url = "http://localhost:9910/BondTools/aligner";
    
    // Inputs for matlab function
    var input = ['record_summary.txt', stack, 'd018'].concat((rc) ? ['radiocarbon'] : []);

    //Use MPS RESTful API to specify params using JSON
    var params = { "nargout": 1,
                   "rhs": input };

    request.open("POST", url);

    //Use MPS RESTful API to set Content-Type
    request.setRequestHeader("Content-Type", "application/json");

    request.onload = function() {   
        //Use MPS RESTful API to check HTTP Status
        if (request.status == 200) {
            // Deserialization: Converting text back into JSON object
            // Response from server is deserialized 
            var result = JSON.parse(request.responseText);
			
            //Use MPS RESTful API to retrieve response in "lhs"
            if('lhs' in result) {  
                document.getElementById("error").innerHTML = "" ;

                // visualizes returned alignment
                visualize(result.lhs[0].mwdata);
            } else { 
                document.getElementById("error").innerHTML = "Error: " + result.error.message; 
            }
        } else { 
            document.getElementById("error").innerHTML = "Error:" + request.statusText; 
        }
    }
    //Serialization: Converting JSON object to text prior to sending request
    request.send(JSON.stringify(params)); 
}
