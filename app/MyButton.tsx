// MyButton.tsx
import React from 'react';

interface MyButtonProps {
  onClick: () => void;
}

const MyButton: React.FC<MyButtonProps> = ({ onClick }) => {
  return (
    <button onClick={onClick}>
      Traducció a Proteïnes
    </button>
  );
};

export default MyButton;