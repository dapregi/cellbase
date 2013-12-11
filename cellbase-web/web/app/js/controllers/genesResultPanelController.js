var genesResult = myApp.controller('genesResult', ['$scope', 'mySharedService', 'CellbaseService', function ($scope, mySharedService, CellbaseService) {

    $scope.toggleTree = []; //array of booleans that will show of hide the elements of the tree
    $scope.genesAndTranscriptsData = {};
    $scope.paginationData = [];
    $scope.biotypes = [];

    $scope.firstGeneId = "";
    $scope.showAll = false;

    $scope.showGenePanel = false;
    $scope.showMoreAndLessGeneData = "+";
    $scope.genePanelMore = false;
    $scope.genePanelStatus = "-";

    $scope.showTranscriptPanel = false;
    $scope.showMoreAndLessTranscriptData = "+";
    $scope.transcriptPanelMore = false;
    $scope.transcriptPanelStatus = "-";

    $scope.showPagination = false;
    $scope.firstPages = false;
    $scope.previousPage = false;
    $scope.nextPage = true;
    $scope.lastPages = true;
    $scope.paginationNumbers = [1, 2, 3];
    $scope.maxNumberPagination;
    $scope.numDataPerPage = 10;
    $scope.showPagination = false;
    $scope.lastPage = 1;
    $scope.disableFirstNumber = true;
    $scope.disableSecondNumber = false;
    $scope.disableThirdNumber = false;

    //========================Pagination==================================
    $scope.goToFirstPage = function () {
        $scope.paginationNumbers[0] = 1;
        $scope.paginationNumbers[1] = 2;
        $scope.paginationNumbers[2] = 3;

        $scope.firstPages = false;
        $scope.previousPage = false;
        $scope.nextPage = true;
        $scope.lastPages = true;

        $scope.disableAndEnablePaginationButtons(1);
        $scope.obtainPaginationLimits(1);
    };
    $scope.goToLastPage = function () {
        $scope.paginationNumbers[0] = $scope.maxNumberPagination - 2;
        $scope.paginationNumbers[1] = $scope.maxNumberPagination - 1;
        $scope.paginationNumbers[2] = $scope.maxNumberPagination;

        $scope.firstPages = true;
        $scope.previousPage = true;
        $scope.nextPage = false;
        $scope.lastPages = false;

        $scope.disableAndEnablePaginationButtons($scope.maxNumberPagination);
        $scope.obtainPaginationLimits($scope.maxNumberPagination);
    };
    $scope.goPreviousPage = function () {
        var page = $scope.lastPage - 1;

        $scope.firstPages = true;
        $scope.previousPage = true;
        $scope.nextPage = true;
        $scope.lastPages = true;

        if (page == 1) {
            $scope.firstPages = false;
            $scope.previousPage = false;

            $scope.paginationNumbers[0] = 1;
            $scope.paginationNumbers[1] = 2;
            $scope.paginationNumbers[2] = 3;
        }
        else if ($scope.paginationNumbers[0] != page && $scope.paginationNumbers[1] != page && $scope.paginationNumbers[2] != page) {
            $scope.paginationNumbers[0] = page - 2;
            $scope.paginationNumbers[1] = page - 1;
            $scope.paginationNumbers[2] = page;
        }
        $scope.disableAndEnablePaginationButtons(page);
        $scope.obtainPaginationLimits(page);
    };
    $scope.goNextPage = function () {
        var page = $scope.lastPage + 1;

        $scope.firstPages = true;
        $scope.previousPage = true;
        $scope.nextPage = true;
        $scope.lastPages = true;

        if (page == $scope.maxNumberPagination) {
            $scope.nextPage = false;
            $scope.lastPages = false;

            $scope.paginationNumbers[0] = page - 2;
            $scope.paginationNumbers[1] = page - 1;
            $scope.paginationNumbers[2] = page;
        }
        else if ($scope.paginationNumbers[0] != page && $scope.paginationNumbers[1] != page && $scope.paginationNumbers[2] != page) {
            $scope.paginationNumbers[0] = page;
            $scope.paginationNumbers[1] = page + 1;
            $scope.paginationNumbers[2] = page + 2;
        }
        $scope.disableAndEnablePaginationButtons(page);
        $scope.obtainPaginationLimits(page);
    };
    $scope.goToNumberPage = function (selectedPage) {
        if (!$scope.simplePagination) {
            if (selectedPage == $scope.maxNumberPagination) {
                $scope.nextPage = false;
                $scope.lastPages = false;
                $scope.firstPages = true;
                $scope.previousPage = true;
            }
            else if (selectedPage == 1) {
                $scope.firstPages = false;
                $scope.previousPage = false;
                $scope.nextPage = true;
                $scope.lastPages = true;
            }
            else {
                $scope.firstPages = true;
                $scope.previousPage = true;
                $scope.nextPage = true;
                $scope.lastPages = true;
            }
        }
        $scope.collapseAllGenesTree();
        $scope.disableAndEnablePaginationButtons(selectedPage);
        $scope.obtainPaginationLimits(selectedPage);
    };
    $scope.disableAndEnablePaginationButtons = function (page) {
        if ($scope.paginationNumbers[0] == page) {
            $scope.disableFirstNumber = true;
            $scope.disableSecondNumber = false;
            $scope.disableThirdNumber = false;
        }
        else if ($scope.paginationNumbers[1] == page) {
            $scope.disableSecondNumber = true;
            $scope.disableFirstNumber = false;
            $scope.disableThirdNumber = false;
        }
        else {
            $scope.disableThirdNumber = true;
            $scope.disableSecondNumber = false;
            $scope.disableFirstNumber = false;
        }
    };
    $scope.obtainPaginationLimits = function (page) {
        $scope.lastPage = page;
        var ini = (page - 1) * $scope.numDataPerPage;
        $scope.paginationData = [];
        var geneId;

        for (var i = ini; i < ini + $scope.numDataPerPage; i++) {
            geneId = Object.keys($scope.genesAndTranscriptsData)[i];
            if (Object.keys($scope.genesAndTranscriptsData)[i] != null) {
                $scope.paginationData.push($scope.genesAndTranscriptsData[geneId]);
            }
        }
    };
    $scope.initPagination = function () {
        $scope.paginationData = [];
        $scope.maxNumberPagination = Math.ceil(Object.keys($scope.genesAndTranscriptsData).length / $scope.numDataPerPage);

        //  0 --> 10
        if (Object.keys($scope.genesAndTranscriptsData).length <= $scope.numDataPerPage) {
            for (var i in $scope.genesAndTranscriptsData) {
                $scope.paginationData.push($scope.genesAndTranscriptsData[i]);
            }
            $scope.showPagination = false;
        }
        // 11 --> 20
        else if (Object.keys($scope.genesAndTranscriptsData).length <= ($scope.numDataPerPage * 2)) {
            $scope.simplePagination = true;

            for (var i = 0; i < $scope.numDataPerPage; i++) {
                geneId = Object.keys($scope.genesAndTranscriptsData)[i];
                if (Object.keys($scope.genesAndTranscriptsData)[i] != null) {
                    $scope.paginationData.push($scope.genesAndTranscriptsData[geneId]);
                }
            }
            $scope.showPagination = true;
            $scope.lastPage = 1;

            $scope.disableFirstNumber = true;
            $scope.disableSecondNumber = false;
            $scope.disableThirdNumber = false;

            $scope.firstPages = false;
            $scope.previousPage = false;
            $scope.nextPage = false;
            $scope.lastPages = false;

            $scope.thirdNumber = false;
            $scope.paginationNumbers = [1, 2];
        }
        // 21 --> ...
        else {
            $scope.simplePagination = false;
            var geneId;

            for (var i = 0; i < $scope.numDataPerPage; i++) {
                geneId = Object.keys($scope.genesAndTranscriptsData)[i];
                if (Object.keys($scope.genesAndTranscriptsData)[i] != null) {
                    $scope.paginationData.push($scope.genesAndTranscriptsData[geneId]);
                }
            }
            $scope.firstPages = false;
            $scope.previousPage = false;
            $scope.nextPage = true;
            $scope.lastPages = true;

            $scope.thirdNumber = true;
            $scope.paginationNumbers = [1, 2, 3];
            $scope.showPagination = true;
            $scope.lastPage = 1;

            $scope.disableFirstNumber = true;
            $scope.disableSecondNumber = false;
            $scope.disableThirdNumber = false;
        }
    };


    $scope.clearAll = function(){
        $scope.showAll = false;
    };
    $scope.clear = function () {
        $scope.showGenePanel = false;
        $scope.showTranscriptPanel = false;
    };
    $scope.setResult = function(){
        $scope.genesFilters = mySharedService.genesIdFilter;
        $scope.biotypeFilters = mySharedService.biotypesFilter;
        $scope.selectedSpecie = mySharedService.selectedSpecies;

        $scope.genesAndTranscriptsData = {};

        var genesIdFilter = [];
        var arrayOfGenes = [];

        //check if there are filters
        if ($scope.biotypeFilters.length != 0) {
            arrayOfGenes = CellbaseService.getGenesAndTranscripts($scope.selectedSpecie.shortName, mySharedService.regionsAndChromosomes, $scope.biotypeFilters);

            for (var i in arrayOfGenes) {
                $scope.genesAndTranscriptsData[arrayOfGenes[i].id] = arrayOfGenes[i];
            }
        }
        if ($scope.genesFilters.length != 0) {
            genesIdFilter = CellbaseService.getGenesAndTranscriptsByIdOrName($scope.selectedSpecie.shortName, $scope.genesFilters);  //obtener los datos

            $scope.checkGeneFilter(genesIdFilter)
        }
        //if there aren't any filters, show all genes data
        if ($scope.biotypeFilters.length == 0 && $scope.genesFilters.length == 0) {
            arrayOfGenes = CellbaseService.getGenesAndTranscripts($scope.selectedSpecie.shortName, mySharedService.regionsAndChromosomes, []);
            //save the data in a hash table
            for (var i in arrayOfGenes) {
                $scope.genesAndTranscriptsData[arrayOfGenes[i].id] = arrayOfGenes[i];
            }
            $scope.getBiotypes();
        }
        $scope.numResults = Object.keys($scope.genesAndTranscriptsData).length;
        $scope.initPagination();
        $scope.clear();
        if($scope.numResults != 0){
            $scope.toggleTree = [];

            for(var i=0;i< 10; i++){
                $scope.toggleTree.push(false);
            }
            $scope.showAll = true;
            $scope.firstGeneId = Object.keys($scope.genesAndTranscriptsData)[0];
            $scope.lastDataShow = Object.keys($scope.genesAndTranscriptsData)[0];
            $scope.selectedGene = CellbaseService.getGenesAllDataById($scope.selectedSpecie.shortName, $scope.lastDataShow);
            //show the informtion of the first gen
            $scope.showSelectedGene(Object.keys($scope.genesAndTranscriptsData)[0], 0);
        }
        else{
            alert("No results with this data");
//            alert("No correct data selected");
            $scope.paginationData = [];
        }
    };
    //save thee correct results and alert the incorrect
    $scope.checkGeneFilter = function(genesIdFilter){
        var genesIdError = [];
        var genesFilters =  $scope.genesFilters.split(",");
        var error = false;

        for(var i in genesIdFilter){
            if(genesIdFilter[i] == undefined){
                genesIdError.push(genesFilters[i]);
                error = true
            }
            else{
                $scope.genesAndTranscriptsData[genesIdFilter[i].id] = (genesIdFilter[i]);
            }
        }
        if(error){
        var messageError = "";
        if(genesIdError.length != 0){
            messageError = genesIdError[0];
            for(var i=1;i<genesIdError.length;i++){
                messageError = messageError + ", " + genesIdError[i];
            }
        }
        messageError = messageError + " incorrect";
        alert(messageError);
        }
    };
    //obtain the list of the biotypes
    $scope.getBiotypes = function () {
        $scope.biotypes = [];
        for (var i in $scope.genesAndTranscriptsData) {
            if ($scope.biotypes.indexOf($scope.genesAndTranscriptsData[i].biotype) == -1) {
                $scope.biotypes.push($scope.genesAndTranscriptsData[i].biotype);
            }
        }
        mySharedService.broadcastGenesBiotypes($scope.biotypes);
    };
    //===================== tree events ========================
    //show gen panel
    $scope.showSelectedGene = function (geneId, index) {
        if($scope.toggleTree[index]){
            $scope.toggleTree[index] = false;
        }
        else{
            $scope.toggleTree[index] = true;
        }
        if ($scope.lastDataShow != geneId) {
            $scope.lastDataShow = geneId;
            $scope.showGenePanel = true;
            $scope.selectedGene = CellbaseService.getGenesAllDataById($scope.selectedSpecie.shortName, geneId);

            $scope.showTranscriptPanel = false;
            $scope.expandAllPanels();
        }
        else {
            if (!$scope.showGenePanel) {
                $scope.showGenePanel = true;
            }
        }
        $scope.selectedTranscripts = $scope.selectedGene.transcripts;
    };
    //show transcripts panel
    $scope.showSelectedTranscript = function (geneId, transcriptName) {
        var transcripts;

        if ($scope.lastDataShow != geneId) {
            $scope.lastDataShow = geneId;
            $scope.showGenePanel = false;
            $scope.selectedGene = CellbaseService.getGenesAllDataById($scope.selectedSpecie.shortName, geneId);
            $scope.expandAllPanels();
        }
        $scope.showTranscriptPanel = true;
        transcripts = $scope.selectedGene.transcripts;
        for (var i in transcripts) {
            if (transcripts[i].name == transcriptName) {
                $scope.selectedTranscript = transcripts[i];
            }
        }
    };

    //show transcripts panel from transcripts table
    $scope.showTanscriptFromTable = function (transcriptName) {
        var transcripts = $scope.selectedGene.transcripts;
        for (var i in transcripts) {
            if (transcripts[i].name == transcriptName) {
                $scope.selectedTranscript = transcripts[i];
            }
        }
        $scope.transcriptInfo = false;
        $scope.showTranscriptPanel = true;
    };

    //Expand/collapse elements in DOM
    $scope.expandAllPanels = function () {
        $scope.geneInfo = false;
        $scope.transcriptInfo = false;
        $scope.genePanelStatus = "-";
        $scope.transcriptPanelStatus = "-";
    };
    $scope.collapseAllPanels = function () {
        $scope.geneInfo = true;
        $scope.transcriptInfo = true;
        $scope.genePanelStatus = "+";
        $scope.transcriptPanelStatus = "+";
    };
    $scope.expandAllGenesTree = function () {
        for(var i in $scope.toggleTree){
            $scope.toggleTree[i] = true;
        }
    };
    $scope.collapseAllGenesTree = function () {
        for(var i in $scope.toggleTree){
            $scope.toggleTree[i] = false;
        }
    };

    //show more info in gen panel
    $scope.showMoreGeneData = function () {
        $scope.genePanelMore = !$scope.genePanelMore;
        if ($scope.showMoreAndLessGeneData == "+") {
            $scope.showMoreAndLessGeneData = "-";
        }
        else {
            $scope.showMoreAndLessGeneData = "+";
        }
    };
    //show more info in transcript panel
    $scope.showMoreTranscriptData = function () {
        $scope.transcriptPanelMore = !$scope.transcriptPanelMore;
        if ($scope.showMoreAndLessTranscriptData == "+") {
            $scope.showMoreAndLessTranscriptData = "-";
        }
        else {
            $scope.showMoreAndLessTranscriptData = "+";
        }
    };

    //show/hide gen panel information
    $scope.openCloseGenePanel = function () {
        if ($scope.genePanelStatus == "+") {
            $scope.genePanelStatus = "-";
        }
        else {
            $scope.genePanelStatus = "+";
        }
    };
    //show/hide transcript panel information
    $scope.openCloseTranscriptPanel = function () {
        if ($scope.transcriptPanelStatus == "+") {
            $scope.transcriptPanelStatus = "-";
        }
        else {
            $scope.transcriptPanelStatus = "+";
        }
    };

    //genesResult div width is the rest of the document
    $scope.getWidth = function () {
        var resultPartWidth = $(document).width() - 220 - 260 - 60;

        console.log(resultPartWidth);
        return  {width: resultPartWidth}
    };
    //tabs
    $scope.goToTab = function () {
        $(function () {
            $('#transcriptsTab a:first').tab('show')
        })
        $('#transcriptsTab a').click(function (e) {
            e.preventDefault()
            $(this).tab('show')
        })
    };

    //--------------EVENTS-------------------
    $scope.$on('clear', function () {
        $scope.clearAll();
    });
    $scope.$on('newSpecie', function () {
        $scope.clearAll();
    });
    $scope.$on('genesClear', function () {
        $scope.clearAll();
    });
    $scope.$on('genesNewResult', function () {
        $scope.setResult();
    });

}]);

genesResult.$inject = ['$scope', 'mySharedService'];