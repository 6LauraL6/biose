// page.tsx
'use client';
import seq from "bionode-seq";
import { useEffect, useState } from "react";
import { select, scaleLinear, scaleBand, scaleOrdinal, schemeCategory10, max, axisBottom, axisLeft } from 'd3';

export default function MainPage() {

  var [dna, setDna] = useState("");
  const [proteinSequence, setProteinSequence] = useState('');
  const [cds_seg, setCDSSeg] = useState('');
  const [currentTime, setCurrentTime] = useState(0);
  const [chartData, setChartData] = useState<number[]>([]);
  const [baseCounts, setBaseCounts] = useState({ A: 0, T: 0, C: 0, G: 0, });
  const [codonTable, setCodonTable] = useState('1');

  // Mostra la hora al carregar la pàgina.
  useEffect(() => {
    fetch('/api/time').then(res => res.json()).then(data => {
      setCurrentTime(data.time);
    });
  }, [])

  // Afegir el càlcul del percentatge GC només si s'ha introduït alguna lletra
  const gcPercentage = dna.length > 0 ? calculateGCPercentage(dna) : 0;
  console.log("Percentatge GC:", gcPercentage);

  let dna_template = seq.reverseComplement(dna);
  let rna = seq.transcribe(dna);

  let seq_prot = ""

  function handleDNA(e: Event) {
    e.preventDefault();
    let input = e.target as HTMLInputElement;
    let value = input.value.toUpperCase(); // converteix a majúscules
    let id = input.id;
    // No apliquis cap filtre a codonTable
    if (id === "codonTable") {
      setCodonTable(value);
      return;
    }
    let filteredValue = value.replace(/[^ACGT]/g, '');
    // console.log("filteredValue="+filteredValue)
    value = seq.reverseComplement(filteredValue);
    setDna(filteredValue);
    // Esborra la lletra d'ADN incorrecta.
    input.value = filteredValue;
  }
  
  function calculateGCPercentage(seq: string): number {
    const gcCount = (seq.match(/[GC]/g) || []).length;
    const totalBases = seq.length;
    return (gcCount / totalBases) * 100;
  }

  // Tradueix la cadena de nucleòtids a aminoacids i fa un recompte de les bases GC.
  const handleTranslateRna = async () => {
    if (!dna || dna.trim() === '') {
      alert("Si us plau, introdueix una seqüència d'ADN.");
      return;
    } else {
      try {
        const response = await fetch('/api/translate', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ mrna_sequence: rna , translation_table : codonTable }),
        });

        if (response.ok) {
          const data = await response.json();
          console.log(data.protein_sequence);
          setProteinSequence(data.protein_sequence);
          setCDSSeg(data.cds_seg);
        } else {
          console.error('Error en la crida POST:', response.statusText);
        }
      } catch (error) {
        console.error('Error en la crida POST:', error);
      }
    }
  };

  // Fa un recompte de les bases amb Javascript i després pinta un gràfic.
  const handleBarPlots = async () => {
    if (!dna || dna.trim() === '') {
      alert("Si us plau, introdueix una seqüència d'ADN abans de generar el gràfic.");
      return;
    }
    const baseCounts = {
      A: (dna.match(/A/g) || []).length,
      T: (dna.match(/T/g) || []).length,
      C: (dna.match(/C/g) || []).length,
      G: (dna.match(/G/g) || []).length,
    };

    const countsArray = Object.values(baseCounts);
    setChartData(countsArray);
  };
  
  // Aquest useEffect actualitza el gràfic cada vegada que es canvia el chartData
  useEffect(() => {
    if (chartData.length > 0) {
      // D3.js code per crear el gràfic de barras
      const svg = select('#chart');
      svg.selectAll('*').remove(); // Esborrar qualsevol gràfic anterior

      const width = 400;
      const height = 300;
      const margin = { top: 20, right: 30, bottom: 40, left: 40 };

      const xScale = scaleBand()
        .domain(Object.keys(baseCounts))
        .range([margin.left, width - margin.right])
        .padding(0.2);

      const yScale = scaleLinear()
        .domain([0, max(chartData) || 1])
        .range([height - margin.bottom, margin.top]);

      const svgChart = svg.append('g')
        .attr('transform', `translate(${margin.left}, 0)`);

      const colorScale = scaleOrdinal(schemeCategory10);

      svgChart.selectAll('rect')
        .data(chartData)
        .enter()
        .append('rect')
        .attr('x', (d, i) => xScale(Object.keys(baseCounts)[i]) || 0)
        .attr('y', d => yScale(d))
        .attr('width', xScale.bandwidth())
        .attr('height', d => height - margin.bottom - yScale(d))
        .attr('fill', (d, i) => colorScale(String(i)));

      svgChart.append('g')
        .attr('transform', `translate(0, ${height - margin.bottom})`)
        .call(axisBottom(xScale));

      svgChart.append('g')
        .attr('transform', `translate(${margin.left}, 0)`)
        .call(axisLeft(yScale));
    }
  }, [chartData, baseCounts]);

  
  return (
    <main className="p-6">
      <h1 className="text-indigo-700 text-xl font-bold">BioActivitat 3 - Expressió gènica</h1>
      <div className="grid grid-cols-2 gap-2 bg-white shadow-md rounded px-8 pt-6 pb-8 mb-4">
        <div className="space-x-2">
          <Input id="dna-coding" pattern="/[^ACGT]/g" readOnly={false} label="DNA Coding Strand" seq={dna} handler={handleDNA} />
          <Input id="codonTable" label="codonTable" pattern="/^[0-9]+$" seq={codonTable} readOnly={false} handler={handleDNA} />
          <Button onClick={handleTranslateRna} buttonText="Tradueix a proteïnes" />
          <Button onClick={handleBarPlots} buttonText="Gràfics de bases" />
        </div>
        <div className="space-x-2">
          <h2 className="text-gray-700 text-xl font-bold mb-2">Resultats.</h2>
          {proteinSequence && <h3 className="block text-gray-700 text-xl pt-5 font-bold">Seqüència de proteïnes: {proteinSequence}</h3>}
          {cds_seg && <h3 className="block text-gray-700 text-xl pt-5 font-bold">Segments: {cds_seg}</h3>}
          {/* Mostra el percentatge GC només si s'ha introduït alguna lletra */}
          {dna.length > 0 && <p>Percentatge GC: {gcPercentage.toFixed(2)}%</p>}
          {/* gràfic amb d3.js */}
          <svg id="chart" width="600" height="400"></svg>
        </div>
      </div>
    </main>
  );
}

// Component camp entrada.

interface InputProps {
  id: string
  label: string
  seq: string
  handler: any
  readOnly: boolean
  pattern: string
}
function Input(props: InputProps) {
  return <div className="mb-4">
    <label className="block text-gray-700 text-sm font-bold mb-2">{props.label}</label>
    <input id={props.id} type="text" readOnly={props.readOnly} pattern={props.pattern} defaultValue={props.seq} 
      onKeyUp={(e) => props.handler(e)} 
      className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline" />
  </div>
}

// Component botó.
// MyButton.tsx
import React from 'react';

interface MyButtonProps {
  onClick: () => void;
  buttonText: string;
}

const Button: React.FC<MyButtonProps> = ({ onClick, buttonText }) => {
  return (
    <button onClick={onClick}>
      {buttonText}
    </button>
  );
};
