// pages/arnm-translation.tsx
import { useState, ChangeEvent } from 'react';

interface TranslationResult {
  translation: string;
}

export function print(data: TranslationResult) {
  console.log(data);
}

const ArnmTranslation: React.FC = () => {
  const [mrnaSequence, setMrnaSequence] = useState<string>('');
  const [translationTable, setTranslationTable] = useState<number>(1);
  const [result, setResult] = useState<string>('');
  // Aquí hi ha una xapussa enorme: http://localhost:5000

  const translate = async () => {
    try {
      const response = await fetch('http://localhost:5000/translate', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          mrna_sequence: mrnaSequence,
          translation_table: translationTable,
        }),
      });

      if (!response.ok) {
        throw new Error('Translation failed');
      }

      const data: TranslationResult = await response.json();
      setResult(JSON.stringify(data, null, 2));
    } catch (error) {
      console.error('Error:', error);
    }
  };

  return (
    <div>
      <h1 style={{ textAlign: 'center' }}>ARNm Translation Form</h1>

      <form>
        <label htmlFor="mrna_sequence">ARNm Sequence:</label>
        <input
          type="text"
          id="mrna_sequence"
          name="mrna_sequence"
          value={mrnaSequence}
          onChange={(e: ChangeEvent<HTMLInputElement>) => setMrnaSequence(e.target.value)}
          required
        />

        <label htmlFor="translation_table">Translation Table:</label>
        <input
          type="number"
          id="translation_table"
          name="translation_table"
          value={translationTable}
          onChange={(e: ChangeEvent<HTMLInputElement>) => setTranslationTable(parseInt(e.target.value, 10))}
          required
        />

        <button type="button" onClick={translate}>
          Translate
        </button>
      </form>

      <pre>{result}</pre>
    </div>
  );
};

export default ArnmTranslation;