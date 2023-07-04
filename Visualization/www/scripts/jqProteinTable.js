function roundVal(value,digits) { return (value>0 || value<0) ? Math.round(value*(10**digits))/(10**digits) : ''; }

function proteinTableCreate() {
  var tableCreated = $.Deferred();

  // load data
  $.getJSON('data/ProteinGeneRecords.json', function(json) {  
    proteinData = json;
    // colNames = Object.keys(proteinData[0]); // ['prtn', 'gene', 'chart', 'desc', 'peptides'] // 2023-06-12
    const colNames = ['prtn', 'gene', 'chart', 'desc','t12_LnLM1', 't12_CFit', 'Npeptides', 'peptides'];
    const colHeaders = {'prtn': 'Protein',  'gene': 'Gene', 'chart': 'Chart', 'desc':'Description', 't12_LnLM1':'t<sub>&frac12; (log-linear)</sub>', 't12_CFit':'t<sub>&frac12; (curve-fit)</sub>', 'Npeptides':'&nbsp;&nbsp;&nbsp;&nbsp;# of<br/>peptides', 'peptides':'Peptides' } ;

    // Create Protein main table
    var table = document.getElementById("tbProtein");

    // build table header
    let tbhead = document.createElement("thead");
    let rowhead = document.createElement("tr");
    colNames.forEach( n => { 
      let th = document.createElement('th'); 
      th.className = colHeaders[n]; 
      th.setAttribute( 'title', "Click to sort table by this column" );
      if ( n=='t12_LnLM1' || n=='t12_CFit' || n=='Npeptides') th.setAttribute("data-sortas", "numeric"); // numeric columns need this
      th.innerHTML = colHeaders[n]; 
      rowhead.appendChild(th); 
    });
    tbhead.appendChild(rowhead);
    table.appendChild(tbhead);

    // build table rows
    let tbbody = document.createElement("tbody");
    proteinData.forEach( proteinrow => {
      const gene = proteinrow['gene']?proteinrow['gene']:proteinrow['prtn'];
      let row = document.createElement("tr");
      row.id = gene;
      tbbody.appendChild(row);

      // add cells
      colNames.forEach( colname => {
        let cell = document.createElement('td');
        cell.className = colname;
        if (colname == "prtn" || colname == "gene") { 
          cell.innerHTML = proteinrow[colname];
          if (proteinrow['chart']>0) {
            cell.style.color = "blue";
            cell.setAttribute( 'title', "Click to load charts" );
            cell.setAttribute( 'onclick', 'showProteinChart(this)');
        }
        } else if (colname == "Npeptides") { 
          cell.innerHTML = proteinrow["peptides"]["Peptide"].length; // need to create hover for peptide list
          cell.style.color = "blue";
          cell.setAttribute( 'title', "Click to toggle show/hide of peptides" );
          cell.setAttribute( 'onclick', 'pepListToggle(this)');
        } else if (colname == "peptides") {
          let div = document.createElement('div');
          div.style.fontSize="small";
          div.style.display="none";
          const thisrow = proteinrow[colname];
          let thisli = [];
          for (let i=0; i<thisrow['Peptide'].length; i++) {
            // let thisval = '';
            let thisval = thisrow['Peptide'][i] + '<br/>LM1: ' + roundVal(thisrow['t12_LnLM1'][i],2)+ ', '+ roundVal(thisrow['r2_LnLM1'][i],4) + '<br/>LM2: ' + roundVal(thisrow['t12_LnLM2'][i],2)+ ', '+ roundVal(thisrow['r2_LnLM2'][i],4) + '<br/>CF: ' + roundVal(thisrow['t12_CFit'][i],2)+ ', '+ roundVal(thisrow['r2_CFit'][i],4) ;
            thisli.push(thisval);
          }
          div.innerHTML = "<ul><li>"+thisli.join("</li><li>") + "</li></ul>";
          div.setAttribute("class", "peptideList");
          cell.appendChild(div);
          cell.setAttribute("class", "peptideList");
          // cell.set
        } else if (colname == "chart") {
          // cell.innerHTML = (proteinrow[colname]>0) ? '&#x2713;' : '&#x274C;'; // 2713/2714
          if (proteinrow[colname]>0) {
            cell.style.color = "green";
            cell.innerHTML = '&#x2713;' ;  // '&#x2713;', '&#x2714;'  
          } else {
            cell.style.color = "red";
            cell.innerHTML = '&#x2717;' ; //  '&#x2717;', '&#x274C;
          }
        } else if (colname == 't12_LnLM1' || colname=='t12_CFit' ) {
          cell.innerHTML = roundVal(proteinrow[colname],2);
        } else { cell.innerHTML = proteinrow[colname]; } // others like 'desc'
        row.appendChild(cell);
      });    
    });
    // dataLoaded.resolve();
    table.appendChild(tbbody);
    let div = $("div#divProteinTable")[0];
    div.innerHTML = null; // remove "Loading table ..."
    div.appendChild(table); // div is an array

    // browser=browserChk(true);
    // if ( !( browser['isChrome'] || browser['isFirefox'] || browser['isIE'] || browser['isEdge'] || browser['isBlink'] )) {
    //   console.log('checked: not chrome, firefox, IE, edge, nor Blink');
    //   // showFinalTable();
    // };

    tableCreated.resolve();
  });

  return tableCreated.promise();
}

function browserChk(log=false) {
  // Safari, for example, complains about Maximum call stack size exceeded with FancyTable iterating thru 10k of data rows. 
  // check broswer: https://stackoverflow.com/questions/9847580/how-to-detect-safari-chrome-ie-firefox-and-opera-browsers
  var isOpera = (!!window.opr && !!opr.addons) || !!window.opera || navigator.userAgent.indexOf(' OPR/') >= 0; // Opera 8.0+
  var isFirefox = typeof InstallTrigger !== 'undefined'; // Firefox 1.0+
  var isSafari = /constructor/i.test(window.HTMLElement) || (function (p) { return p.toString() === "[object SafariRemoteNotification]"; })(!window['safari'] || (typeof safari !== 'undefined' && window['safari'].pushNotification)); // Safari 3.0+ "[object HTMLElementConstructor]" 
  var isIE = /*@cc_on!@*/false || !!document.documentMode; // Internet Explorer 6-11
  var isEdge = !isIE && !!window.StyleMedia; // Edge 20+
  var isChrome = !!window.chrome; // !!window.chrome && (!!window.chrome.webstore || !!window.chrome.runtime); // Chrome 1 - 79
  var isEdgeChromium = isChrome && (navigator.userAgent.indexOf("Edg") != -1); // Edge (based on chromium) detection
  var isBlink = (isChrome || isOpera) && !!window.CSS;  // Blink engine detection

  if (log) { console.log(' isOpera: '+isOpera+'\n isFirefox: '+isFirefox+'\n isSafari: '+isSafari+'\n isIE: '+isIE+'\n isEdge: '+isEdge+'\n isChrome: '+isChrome+'\n isEdgeChromium: '+isEdgeChromium+'\n isBlink: '+isBlink); }

  return {isOpera:isOpera, isFirefox:isFirefox, isSafari:isSafari, isIE:isIE, isEdge:isEdge, isChrome:isChrome, isEdgeChromium:isEdgeChromium, isBlink:isBlink}
}

function showProteinChart(elm) {
  var row = elm.parentNode;
  const gene = row.id; // event.target.attributes['proteinlink'].value;

  // set focus 
  const position =  $("div#proteinChart iframe").position();
  scroll(0,position.top);

  var proteinframe = $("div#proteinChart iframe")[0];
  var peptideframe = $("div#peptideChart iframe")[0];
  const proteinPage = "/charts/proteinLevel/RelAbundance_Gene-"+gene+".html";
  const peptidePage = "/charts/peptideLevel/RelAbundance_Gene-"+gene+"-peptides.html";

  proteinframe.setAttribute('src', proteinPage);
  peptideframe.setAttribute('src', peptidePage);

  // change the pop-up links
  var proteinpopup = $("div#proteinChart a:first")[0];
  var peptidepopup = $("div#peptideChart a:first")[0];
  proteinpopup.setAttribute('href', proteinPage);
  peptidepopup.setAttribute('href', peptidePage);

  return null;
}

function applyFancyTable() {
  // var fancyTableCreated = $.Deferred();
  proteinTable = $("table#tbProtein").fancyTable({
    sortColumn:1,
    pagination: true,
    perPage: 50,
    globalSearch: true,
    globalSearchExcludeColumns: [3,5,6,7]
  });
  // fancyTableCreated.resolve();
  return;
  // return fancyTableCreated.promise();
}

function pepListToggle(elm) {
  var row = elm.parentNode;
  var tdPepList = row.querySelector("div.peptideList");
  var curStatus = tdPepList.style.display;
  if (curStatus==="none") { 
    $("td > div.peptideList").hide();
    tdPepList.style.display = "block"; 
  } else {
    tdPepList.style.display = "none"; 
  }
  return null;
}

function showFinalTable() {
  browser=browserChk(true);
  if ( !( browser['isChrome'] || browser['isFirefox'] || browser['isIE'] || browser['isEdge'] || browser['isBlink'] )) {
    console.log('checked: not chrome, firefox, IE, edge, nor Blink');
    return null;
  };
  
  let t1 = document.getElementById('loadingMsg');
  t1.innerHTML = "&nbsp;&nbsp;Sortable columns (Click on a column header to sort)";
  let t2 = document.getElementById('divProteinTable');
  t2.style.display = 'block';
  return null;
}

